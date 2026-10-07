import test from 'node:test';
import assert from 'node:assert/strict';
import { resolveCompound } from '../dist/chemistry.js';
import { renderCondensed2D } from '../dist/condensed2d.js';

const attributes = tag => Object.fromEntries([...tag.matchAll(/([\w-]+)="([^"]*)"/g)].map(match => [match[1], match[2]]));
const texts = svg => [...svg.matchAll(/<text\b([^>]+)>(.*?)<\/text>/gs)].map(match => ({ ...attributes(match[1]), markup: match[2] }));
const lines = svg => [...svg.matchAll(/<line\b([^>]+)\/>/g)].map(match => attributes(match[1]));
const groups = svg => [...svg.matchAll(/<g class="condensed2d-branch"([^>]+)>/g)].map(match => attributes(match[1]));
const mainLabels = svg => texts(svg).filter(text => text.class.includes('main-chain-label'));
const branchLabels = svg => texts(svg).filter(text => text.class.includes('branch-label'));
const bounds = svg => attributes(svg.match(/^<svg\b([^>]+)>/)[1]);
const deepFreeze = value => {
  if (value && typeof value === 'object') { Object.values(value).forEach(deepFreeze); Object.freeze(value); }
  return value;
};

test('requested octane has all eight explicit groups with ethyl at 3 and tert-butyl at 5', async () => {
  const compound = await resolveCompound('5-terz-butil-3-etilottano');
  const svg = renderCondensed2D(compound, { width: 320, height: 240, numbering: true });
  const chain = mainLabels(svg), branches = groups(svg);
  assert.equal(compound.formula, 'C14H30');
  assert.deepEqual(chain.map(atom => atom['data-label']), ['CH3', 'CH2', 'CH', 'CH2', 'CH', 'CH2', 'CH2', 'CH3']);
  assert.deepEqual(branches.map(branch => [+branch['data-parent-locant'], branch['data-label']]), [[3, 'CH2-CH3'], [5, 'C(CH3)3']]);
  assert.equal(branches[0]['data-parent'], String(compound.backbone[2]));
  assert.equal(branches[1]['data-parent'], String(compound.backbone[4]));
  assert.deepEqual(texts(svg).filter(atom => atom.class === 'atom-number').map(atom => +atom['data-locant']), [1, 2, 3, 4, 5, 6, 7, 8]);
  assert.equal(chain.filter(atom => atom['data-label'] === 'CH2').length, 4);
  assert.ok(branchLabels(svg).every(branch => +branch.y < +chain[0].y));
  assert.match(svg, /baseline-shift="sub">3<\/tspan>/);
  assert.equal(bounds(svg).height, '240');
  assert.ok(+bounds(svg).width > 320, 'The chain retains readable intrinsic width on a phone.');
  // Each rendered compact branch accounts for its actual subtree, without losing atoms or hydrogens.
  const represented = [...chain.map(atom => +atom['data-atom']), ...branches.flatMap(branch => branch['data-atoms'].split(',').map(Number))];
  assert.equal(new Set(represented).size, compound.atoms.length);
  assert.equal(represented.length, compound.atoms.length);
  assert.equal(represented.reduce((total, index) => total + compound.atoms[index].hydrogens, 0), 30);
});

test('branch bonds meet the first CH2 or C group, not the middle of a whole substituent label', async () => {
  const compound = await resolveCompound('5-terz-butil-3-etilottano');
  const svg = renderCondensed2D(compound, { width: 760, height: 240 });
  const branches = groups(svg), labels = branchLabels(svg), bonds = lines(svg);
  for (const branch of branches) {
    const bond = bonds.find(line => line['data-bond'] === `${branch['data-parent']}-${branch['data-root']}`);
    const label = labels.find(text => text['data-atom'] === branch['data-root']);
    assert.equal(bond.x1, bond.x2, 'Each normal substituent is vertically attached.');
    // CH2 has a 22px C and H plus a 14px subscript; C is one 22px glyph.
    const rootHalf = branch['data-label'].startsWith('CH2') ? 17.4 : 6.6;
    assert.ok(Math.abs(+label.x + rootHalf - +bond.x2) < .001);
    assert.ok(+bond.y2 > +label.y && +bond.y1 < +mainLabels(svg)[0].y);
    assert.ok(+bond.y1 > +bond.y2, 'A branch bond has positive visible length.');
  }
});

test('the initial decane example retains an isobutyl rooted CH2-CH(CH3)2 group', async () => {
  const compound = await resolveCompound('5-isobutil decano');
  const svg = renderCondensed2D(compound, { width: 320, height: 240, numbering: true });
  assert.equal(mainLabels(svg).length, 10);
  assert.deepEqual(groups(svg).map(branch => [+branch['data-parent-locant'], branch['data-label']]), [[5, 'CH2-CH(CH3)2']]);
  assert.equal(groups(svg)[0]['data-parent'], String(compound.backbone[4]));
  assert.equal(groups(svg)[0]['data-atoms'].split(',').length, 4);
  assert.equal(mainLabels(svg).filter(atom => atom['data-label'] === 'CH2').length, 7);
  assert.ok(!svg.includes('(CH2)3') && !svg.includes('(CH2)4'));
});

test('two substituents on one carbon occupy opposite sides and preserve quaternary C', async () => {
  const compound = await resolveCompound('2,2-dimetilbutano');
  const svg = renderCondensed2D(compound, { width: 320, height: 240, numbering: true });
  assert.deepEqual(mainLabels(svg).map(atom => atom['data-label']), ['CH3', 'C', 'CH2', 'CH3']);
  assert.deepEqual(groups(svg).map(branch => [+branch['data-parent-locant'], branch['data-label'], branch['data-side']]), [[2, 'CH3', 'above'], [2, 'CH3', 'below']]);
  const [above, below] = branchLabels(svg), mainY = +mainLabels(svg)[1].y;
  assert.ok(+above.y < mainY && +below.y > mainY);
  const mainCenter = +mainLabels(svg)[1].x + 6.6;
  const number = texts(svg).find(atom => atom.class === 'atom-number' && atom['data-locant'] === '2');
  assert.notEqual(+number.x, mainCenter, 'Numbering stays clear of the lower branch bond.');
});

test('double and triple backbone bonds retain line multiplicity and hydrogen counts', async () => {
  const alkene = renderCondensed2D(await resolveCompound('but-2-ene'), { width: 320, height: 160 });
  assert.deepEqual(mainLabels(alkene).map(atom => atom['data-label']), ['CH3', 'CH', 'CH', 'CH3']);
  assert.equal(lines(alkene).filter(line => line['data-order'] === '2').length, 2);
  const alkyne = renderCondensed2D(await resolveCompound('pent-1-ino'), { width: 320, height: 160 });
  assert.deepEqual(mainLabels(alkyne).map(atom => atom['data-label']), ['CH', 'C', 'CH2', 'CH2', 'CH3']);
  assert.equal(lines(alkyne).filter(line => line['data-order'] === '3').length, 3);
});

test('heteroatom substituents preserve their multiple attachment bond', () => {
  const compound = { name: 'Carbonile', atoms: [{ element: 'C', hydrogens: 3 }, { element: 'C', hydrogens: 0 }, { element: 'C', hydrogens: 3 }, { element: 'O', hydrogens: 0 }], bonds: [{ from: 0, to: 1, order: 1 }, { from: 1, to: 2, order: 1 }, { from: 1, to: 3, order: 2 }], backbone: [0, 1, 2] };
  const svg = renderCondensed2D(compound, { width: 320, height: 240 });
  assert.deepEqual(mainLabels(svg).map(atom => atom['data-label']), ['CH3', 'C', 'CH3']);
  assert.equal(groups(svg)[0]['data-label'], 'O');
  assert.equal(lines(svg).filter(line => line['data-bond'] === '1-3').length, 2);
});

test('rings and unidentified parent chains use the complete semideveloped graph', async () => {
  for (const name of ['cicloesano', 'caffeina']) {
    const compound = await resolveCompound(name), svg = renderCondensed2D(compound);
    assert.equal(bounds(svg)['data-fallback'], 'true', name);
    assert.equal(bounds(svg)['data-layout'], 'graph', name);
    assert.match(bounds(svg)['aria-label'], /semisviluppata del grafo/);
    assert.match(svg, /anelli e tutti i loro legami restano visibili/);
    assert.equal(new Set(lines(svg).map(line => line['data-bond'])).size, compound.bonds.length, name);
    assert.equal(texts(svg).filter(atom => atom.class === 'atom-label').length, compound.atoms.length, name);
    assert.equal(mainLabels(svg).length, 0);
  }
  const invalid = await resolveCompound('pentano');
  invalid.backbone = [0, 2, 4];
  assert.equal(bounds(renderCondensed2D(invalid))['data-fallback'], 'true', 'Disconnected backbone indices are never fabricated into bonds.');
});

test('names and unusual element text are escaped in SVG and all graph data stays immutable', async () => {
  const compound = await resolveCompound('5-terz-butil-3-etilottano');
  compound.name = '"<script>alert(1)</script>&';
  compound.atoms[compound.backbone[0]].element = '<img onerror="alert(1)">&';
  const before = JSON.stringify(compound);
  deepFreeze(compound);
  const svg = renderCondensed2D(compound, { width: 320, height: 240, numbering: true });
  assert.ok(!svg.includes('<script>') && !svg.includes('<img '));
  assert.ok(svg.includes('&lt;script&gt;') && svg.includes('&lt;img onerror=&quot;alert(1)&quot;&gt;&amp;'));
  assert.equal(JSON.stringify(compound), before);
  const caffeine = deepFreeze(await resolveCompound('caffeina'));
  assert.doesNotThrow(() => renderCondensed2D(caffeine));
});

test('long chains maintain forward visible bonds and coordinates within intrinsic bounds', async () => {
  const compound = await resolveCompound('eicosano');
  const svg = renderCondensed2D(compound, { width: 320, height: 120, numbering: true });
  assert.equal(mainLabels(svg).length, 20);
  assert.equal(mainLabels(svg).filter(atom => atom['data-label'] === 'CH2').length, 18);
  assert.ok(+bounds(svg).width > 1500);
  for (const line of lines(svg)) assert.ok(+line.x2 > +line.x1, 'Every horizontal bond has forward positive length.');
  const width = +bounds(svg).width, height = +bounds(svg).height;
  for (const item of [...texts(svg), ...lines(svg)]) {
    for (const key of ['x', 'x1', 'x2']) if (key in item) assert.ok(+item[key] >= 0 && +item[key] <= width, `${key} within viewBox`);
    for (const key of ['y', 'y1', 'y2']) if (key in item) assert.ok(+item[key] >= 0 && +item[key] <= height, `${key} within viewBox`);
  }
});
