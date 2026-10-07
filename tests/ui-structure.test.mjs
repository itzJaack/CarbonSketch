import test from 'node:test';
import assert from 'node:assert/strict';
import { resolveCompound } from '../dist/chemistry.js';
import { renderTopology, renderLewis } from '../dist/renderers.js';

const atomLabels = svg => [...svg.matchAll(/<text\b[^>]*data-element="([^"]+)"/g)].map(match => match[1]);
const uniqueBonds = svg => new Set([...svg.matchAll(/data-bond="([^"]+)"/g)].map(match => match[1]));
const lonePairs = svg => [...svg.matchAll(/class="lone-pair"/g)].length;

test('the requested Lewis diagram and exported SVG retain all 44 atoms and 43 bonds', async () => {
  const compound = await resolveCompound('5-isobutil decano');
  for (const options of [{}, { width: 1500, height: 850 }]) {
    const svg = renderLewis(compound, options);
    const labels = atomLabels(svg);
    assert.equal(labels.filter(element => element === 'C').length, 14);
    assert.equal(labels.filter(element => element === 'H').length, 30);
    assert.equal(uniqueBonds(svg).size, 43);
    assert.equal(lonePairs(svg), 0);
    assert.doesNotMatch(svg, /NaN|Infinity|undefined/);
    assert.match(svg, /xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
  }
});

test('topological numbering follows the ten-carbon parent without hiding branch bonds', async () => {
  const compound = await resolveCompound('5-isobutil decano');
  const svg = renderTopology(compound, { numbering: true });
  assert.equal(uniqueBonds(svg).size, 13);
  const numbers = [...svg.matchAll(/class="atom-number"[^>]*>(\d+)<\/text>/g)].map(match => Number(match[1]));
  assert.deepEqual(numbers, Array.from({ length: 10 }, (_, index) => index + 1));
  assert.doesNotMatch(svg, /class="atom-label"/);
});

test('Lewis diagrams show free electron pairs separately from explicit hydrogen atoms', async () => {
  const cases = [
    ['acqua', 3, 2, 2],
    ['etanolo', 9, 8, 2],
    ['acido acetico', 8, 7, 4],
  ];
  for (const [name, atoms, bonds, pairs] of cases) {
    const svg = renderLewis(await resolveCompound(name));
    assert.equal(atomLabels(svg).length, atoms, name);
    assert.equal(uniqueBonds(svg).size, bonds, name);
    assert.equal(lonePairs(svg), pairs, name);
    assert.equal([...svg.matchAll(/class="electron"/g)].length, pairs * 2, name);
  }
});

test('a benzene Kekule drawing keeps all six carbons, hydrogens and alternating double bonds', async () => {
  const compound = await resolveCompound('benzene');
  const svg = renderLewis(compound);
  assert.equal(atomLabels(svg).filter(element => element === 'C').length, 6);
  assert.equal(atomLabels(svg).filter(element => element === 'H').length, 6);
  assert.equal(uniqueBonds(svg).size, 12);
  assert.equal([...svg.matchAll(/data-order="2"/g)].length, 6);
  assert.equal(lonePairs(svg), 0);
});

test('formula exports escape compound names and remain standalone SVG documents', async () => {
  const compound = await resolveCompound('metano');
  compound.name = '<script>alert("x")</script> & molecule';
  for (const render of [renderTopology, renderLewis]) {
    const svg = render(compound);
    assert.doesNotMatch(svg, /<script>/);
    assert.match(svg, /&lt;script&gt;/);
    assert.match(svg, /<title>/);
    assert.match(svg, /<desc>/);
    assert.match(svg, /<style>/);
    assert.ok(svg.endsWith('</svg>'));
  }
});
