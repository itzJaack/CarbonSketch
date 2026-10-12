import test from 'node:test';
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { resolveCompound, compoundFromGraph, buildLewisGraph } from '../dist/chemistry.js';
import { parseCondensedFormula } from '../dist/condensed-editor.js';
import { renderTopology, renderLewis, renderCondensed } from '../dist/renderers.js';
import { renderCondensed2D } from '../dist/condensed2d.js';

function treeFingerprint(compound) {
  assert.equal(compound.bonds.length, compound.atoms.length - 1);
  const graph = compound.atoms.map(() => []);
  compound.bonds.forEach(bond => {
    graph[bond.from].push({ next: bond.to, order: bond.order });
    graph[bond.to].push({ next: bond.from, order: bond.order });
  });
  function rooted(index, parent = -1) {
    const atom = compound.atoms[index];
    const children = graph[index].filter(bond => bond.next !== parent)
      .map(bond => `${bond.order}:${rooted(bond.next, index)}`).sort();
    return `${atom.element}[${atom.hydrogens || 0},${atom.charge || 0}](${children.join(',')})`;
  }
  return compound.atoms.map((_, index) => rooted(index)).sort()[0];
}

function longestCarbonChain(compound) {
  const graph = compound.atoms.map(() => []);
  compound.bonds.forEach(bond => { graph[bond.from].push(bond.to); graph[bond.to].push(bond.from); });
  let longest = 0;
  function visit(index, parent, length) {
    longest = Math.max(longest, length);
    graph[index].filter(next => next !== parent).forEach(next => visit(next, index, length + 1));
  }
  compound.atoms.forEach((_, index) => visit(index, -1, 1));
  return longest;
}

test('names round-trip through the text editor without changing their molecular connectivity', async () => {
  for (const name of [
    '5-isobutil decano',
    '5-terz-butil-3-etilottano',
    '4,5dietil2,3,4,5tetraetilnonano',
    '4-sec-butilottano',
    '2,2-dimetilbutano',
    '3-metiles-1-ene',
    'etanolo',
    'acetone',
    'acido acetico',
  ]) {
    const named = await resolveCompound(name);
    const drawn = compoundFromGraph(parseCondensedFormula(named.condensed));
    assert.equal(drawn.formula, named.formula, name);
    assert.ok(Math.abs(drawn.mass - named.mass) < 0.01, name);
    assert.equal(treeFingerprint(drawn), treeFingerprint(named), name);
    assert.equal(buildLewisGraph(drawn).atoms.length, buildLewisGraph(named).atoms.length, name);
  }
});

test('the six-ethyl nonane input conserves its graph and does not claim a preferred IUPAC name', async () => {
  const compound = await resolveCompound('4,5dietil2,3,4,5tetraetilnonano');
  assert.equal(compound.formula, 'C21H44');
  assert.equal(compound.atoms.length, 21);
  assert.equal(compound.backbone.length, 9);
  assert.equal(longestCarbonChain(compound), 10);
  assert.equal(compound.canonicalNameVerified, false);
  assert.match(compound.normalizedName, /2,3,4,4,5,5-/);
  assert.ok(compound.notes.some(note => /non è la più lunga/.test(note)));
});

test('a written molecular group is checked again after an editor changes its bond order', async () => {
  const graph = parseCondensedFormula('CH3-CH2-OH');
  graph.bonds[0].order = 2;
  assert.throws(() => compoundFromGraph(graph), /idrogen|valenz|legam/i);
  assert.equal((await resolveCompound('etanolo')).formula, 'C2H6O');
});

test('all structural drawings and explicit Lewis data agree for an editor-created branch', () => {
  const compound = compoundFromGraph(parseCondensedFormula('CH3-CH2-CH(CH2-CH3)-CH2-CH(C(CH3)3)-CH2-CH2-CH3'));
  assert.equal(compound.formula, 'C14H30');
  assert.equal(compound.atoms.length, 14);
  const explicit = buildLewisGraph(compound);
  assert.equal(explicit.atoms.length, 44);
  assert.equal(explicit.bonds.length, 43);
  for (const render of [renderTopology, renderLewis, renderCondensed, renderCondensed2D]) {
    const svg = render(compound);
    assert.match(svg, /^<svg\b/);
    assert.doesNotMatch(svg, /NaN|Infinity|undefined/);
    assert.ok(svg.endsWith('</svg>'));
  }
  const svg = renderCondensed2D(compound);
  const atoms = [...svg.matchAll(/class="atom-label main-chain-label"[^>]*data-atom="(\d+)"/g)].map(match => Number(match[1]));
  const branches = [...svg.matchAll(/data-atoms="([\d,]+)"/g)].flatMap(match => match[1].split(',').map(Number));
  assert.equal(new Set([...atoms, ...branches]).size, 14);
  assert.equal(atoms.length + branches.length, 14);
});

test('project Pages URLs resolve the complete local asset and module graph beneath the project directory', async () => {
  const base = new URL('https://itzjaack.github.io/CarbonSketch/');
  const entry = await readFile(new URL('../dist/index.html', import.meta.url), 'utf8');
  const dependencies = [...entry.matchAll(/(?:src|href)="([^"#][^"]*)"/g)].map(match => match[1]);
  for (const dependency of dependencies) {
    if (/^https?:/.test(dependency)) continue;
    const resolved = new URL(dependency, base);
    assert.ok(resolved.pathname.startsWith('/CarbonSketch/'), dependency);
    if (dependency !== './') await readFile(new URL(`../dist/${dependency.replace(/^\.\//, '')}`, import.meta.url));
  }
  const modules = ['app.js', 'chemistry.js', 'renderers.js', 'condensed2d.js', 'condensed-editor.js'];
  for (const module of modules) {
    const content = await readFile(new URL(`../dist/${module}`, import.meta.url), 'utf8');
    for (const match of content.matchAll(/from\s+['"]([^'"]+)['"]/g)) {
      const resolved = new URL(match[1], new URL(module, base));
      assert.ok(resolved.pathname.startsWith('/CarbonSketch/'), `${module}: ${match[1]}`);
      await readFile(new URL(`../dist/${resolved.pathname.slice('/CarbonSketch/'.length)}`, import.meta.url));
    }
  }
});
