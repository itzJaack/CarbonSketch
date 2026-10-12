import test from 'node:test';
import assert from 'node:assert/strict';
import { parseCondensedFormula } from '../dist/condensed-editor.js';
import { compoundFromGraph } from '../dist/chemistry.js';

function neighbors(graph, index) {
  return graph.bonds.flatMap(bond => bond.from === index ? [{ index: bond.to, order: bond.order }] : bond.to === index ? [{ index: bond.from, order: bond.order }] : []);
}
function hydrogens(graph) { return graph.atoms.reduce((total, atom) => total + atom.hydrogens, 0); }

test('the requested drawn octane preserves eight main groups, two branch locants and C14H30', () => {
  const graph = parseCondensedFormula('CH3-CH2-CH(CH2-CH3)-CH2-CH(C(CH3)3)-CH2-CH2-CH3');
  const compound = compoundFromGraph(graph), main = new Set(graph.backbone);
  assert.equal(compound.formula, 'C14H30');
  assert.equal(graph.atoms.length, 14);
  assert.equal(graph.bonds.length, 13);
  assert.equal(graph.backbone.length, 8);
  assert.deepEqual(graph.backbone.map(index => graph.atoms[index].hydrogens), [3, 2, 1, 2, 1, 2, 2, 3]);
  const ethyl = neighbors(graph, graph.backbone[2]).filter(next => !main.has(next.index));
  assert.equal(ethyl.length, 1);
  assert.equal(graph.atoms[ethyl[0].index].hydrogens, 2);
  const terminal = neighbors(graph, ethyl[0].index).filter(next => next.index !== graph.backbone[2]);
  assert.equal(terminal.length, 1);
  assert.equal(graph.atoms[terminal[0].index].hydrogens, 3);
  const tert = neighbors(graph, graph.backbone[4]).filter(next => !main.has(next.index));
  assert.equal(tert.length, 1);
  assert.equal(graph.atoms[tert[0].index].hydrogens, 0);
  const methyls = neighbors(graph, tert[0].index).filter(next => next.index !== graph.backbone[4]);
  assert.equal(methyls.length, 3);
  methyls.forEach(next => assert.equal(graph.atoms[next.index].hydrogens, 3));
  assert.equal(hydrogens(graph), 30);
});

test('nested isobutyl formula and repeated methylene segments conserve the original graph', () => {
  const graph = parseCondensedFormula('CH3-(CH2)3-CH(CH2-CH(CH3)2)-(CH2)4-CH3');
  assert.equal(compoundFromGraph(graph).formula, 'C14H30');
  assert.equal(graph.backbone.length, 10);
  assert.equal(graph.atoms[graph.backbone[4]].hydrogens, 1);
  const main = new Set(graph.backbone), root = neighbors(graph, graph.backbone[4]).find(next => !main.has(next.index)).index;
  assert.equal(graph.atoms[root].hydrogens, 2);
  const fork = neighbors(graph, root).find(next => next.index !== graph.backbone[4]).index;
  assert.equal(graph.atoms[fork].hydrogens, 1);
  assert.equal(neighbors(graph, fork).filter(next => next.index !== root).length, 2);
  graph.atoms.forEach(atom => ['x', 'y', 'z', 'x2d', 'y2d'].forEach(key => assert.ok(Number.isFinite(atom[key]))));
});

test('terminal repeated CH3 branches produce three distinct methyl atoms on a quaternary carbon', () => {
  const graph = parseCondensedFormula('CH3-C(CH3)3');
  assert.equal(compoundFromGraph(graph).formula, 'C5H12');
  assert.equal(graph.atoms[1].hydrogens, 0);
  assert.equal(neighbors(graph, 1).length, 4);
  assert.equal(graph.atoms.filter(atom => atom.hydrogens === 3).length, 4);
});

test('implicit separators, ordinary spaces and Unicode subscripts/dashes represent the same structure', () => {
  for (const formula of ['CH3CH2CH3', 'CH3-CH2-CH3', ' CH₃ – CH₂ – CH₃ ', 'CH3(CH2)1CH3']) {
    const graph = parseCondensedFormula(formula);
    assert.equal(compoundFromGraph(graph).formula, 'C3H8', formula);
    assert.equal(graph.atoms.length, 3, formula);
  }
});

test('parenthesized repeated CH2 is a chain even without a separating dash', () => {
  const graph = parseCondensedFormula('CH3(CH2)4CH3');
  assert.equal(compoundFromGraph(graph).formula, 'C6H14');
  assert.equal(graph.backbone.length, 6);
  assert.deepEqual(graph.backbone.map(index => neighbors(graph, index).length), [1, 2, 2, 2, 2, 1]);
});

test('unsaturated groups preserve double, triple and carbonyl bonds with exact H counts', () => {
  for (const [formula, expected, order] of [['CH3-CH=CH-CH3', 'C4H8', 2], ['CH≡C-CH2-CH3', 'C4H6', 3], ['CH3-C(=O)-CH3', 'C3H6O', 2], ['CH3-CH2-OH', 'C2H6O', 1], ['CH3-C#N', 'C2H3N', 3]]) {
    const graph = parseCondensedFormula(formula);
    assert.equal(compoundFromGraph(graph).formula, expected, formula);
    assert.ok(graph.bonds.some(bond => bond.order === order), formula);
  }
});

test('written invalid H counts are rejected instead of being repaired to a different molecule', () => {
  for (const formula of ['CH2-CH2', 'CH3-CH3-CH3', 'CH3=CH3', 'CH3-C(CH3)4', 'CH3-O-CH2', 'CH5', 'CH0', 'CH3-(CH3)2-CH3']) {
    assert.throws(() => parseCondensedFormula(formula), Error, formula);
  }
});

test('empty input, unmatched parentheses and incomplete bonds produce actionable errors', () => {
  for (const formula of ['', ' ', null, 'CH3-', 'CH3--CH3', 'CH3-()', 'CH3-CH(CH3', 'CH3-CH3)', 'CH3-CH3.', '<script>CH4</script>', 'C6H6']) {
    assert.throws(() => parseCondensedFormula(formula), Error, String(formula));
  }
});

test('repetition limits bound the total graph and fail before building huge input', () => {
  for (const formula of ['CH3-(CH2)0-CH3', 'CH3-(CH2)101-CH3', 'CH3-(CH2)99-CH3', 'CH3-(CH2)99999999999999999-CH3', '('.repeat(22) + 'CH4' + ')'.repeat(22)]) {
    assert.throws(() => parseCondensedFormula(formula), Error, formula.slice(0, 50));
  }
  assert.equal(parseCondensedFormula('CH3-(CH2)98-CH3').atoms.length, 100);
});

test('the shared drawing validator rejects duplicate/self bonds and preserves H errors after editing', () => {
  const original = parseCondensedFormula('CH3-CH2-CH3');
  const duplicate = structuredClone(original); duplicate.bonds.push({ ...duplicate.bonds[0] });
  assert.throws(() => compoundFromGraph(duplicate), Error);
  const self = structuredClone(original); self.bonds.push({ from: 0, to: 0, order: 1 });
  assert.throws(() => compoundFromGraph(self), Error);
  const disconnected = structuredClone(original); disconnected.bonds.pop();
  assert.throws(() => compoundFromGraph(disconnected), Error);
  const edited = structuredClone(original); edited.bonds[0].order = 2;
  assert.throws(() => compoundFromGraph(edited), Error, 'Changing a bond requires changing written H groups too.');
  assert.deepEqual(original.atoms.map(atom => atom.hydrogens), [3, 2, 3]);
});

test('input graph is not mutated while validating a parsed formula', () => {
  const graph = parseCondensedFormula('CH3-CH(CH3)-CH3'), before = JSON.stringify(graph);
  const freeze = value => { if (value && typeof value === 'object') { Object.values(value).forEach(freeze); Object.freeze(value); } return value; };
  assert.equal(compoundFromGraph(freeze(graph)).formula, 'C4H10');
  assert.equal(JSON.stringify(graph), before);
});
