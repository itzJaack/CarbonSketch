import test from 'node:test';
import assert from 'node:assert/strict';
import { resolveCompound, buildLewisGraph } from '../dist/chemistry.js';
import { renderTopology, renderLewis, renderCondensed, create3DViewer } from '../dist/renderers.js';

const attributes = tag => Object.fromEntries([...tag.matchAll(/([\w-]+)="([^"]*)"/g)].map(match => [match[1], match[2]]));
const lines = svg => [...svg.matchAll(/<line\b([^>]+)\/>/g)].map(match => attributes(match[1]));
const labels = svg => [...svg.matchAll(/<text class="atom-label"([^>]+)>(.*?)<\/text>/gs)].map(match => ({ ...attributes(match[1]), label: match[2] }));

function assertBounds(svg) {
  const [, , , width, height] = svg.match(/viewBox="(\S+) (\S+) (\S+) (\S+)"/);
  for (const match of svg.matchAll(/<(?:line|text|circle)\b([^>]+)>?/g)) {
    const attrs = attributes(match[1]);
    for (const key of ['x', 'x1', 'x2', 'cx']) if (key in attrs) assert.ok(+attrs[key] >= 0 && +attrs[key] <= +width, `${key} outside SVG bounds`);
    for (const key of ['y', 'y1', 'y2', 'cy']) if (key in attrs) assert.ok(+attrs[key] >= 0 && +attrs[key] <= +height, `${key} outside SVG bounds`);
  }
}

test('requested compound has all 14 C, 30 H and 43 bonds in Lewis; only 10 parent-chain locants', async () => {
  const compound = await resolveCompound('5-isobutil decano');
  const lewis = renderLewis(compound), topology = renderTopology(compound, { numbering: true });
  const atoms = labels(lewis);
  assert.equal(atoms.filter(atom => atom['data-element'] === 'C').length, 14);
  assert.equal(atoms.filter(atom => atom['data-element'] === 'H').length, 30);
  assert.equal(lines(lewis).length, 43);
  assert.equal(lines(topology).length, 13);
  const numbers = [...topology.matchAll(/<text class="atom-number"([^>]+)>(\d+)<\/text>/g)].map(match => ({ ...attributes(match[1]), number: +match[2] }));
  assert.deepEqual(numbers.map(atom => atom.number), [1,2,3,4,5,6,7,8,9,10]);
  assert.equal(numbers.find(atom => atom.number === 5)['data-atom'], '0');
  assertBounds(lewis); assertBounds(topology);
});

test('Lewis displays correct oxygen/nitrogen lone pairs and retains multiple bonds', async () => {
  for (const [name, atomCount, lonePairs] of [['etanolo',9,2],['acetone',10,2],['acido acetico',8,4],['caffeina',24,8],['acqua',3,2]]) {
    const compound = await resolveCompound(name), svg = renderLewis(compound);
    assert.equal(labels(svg).length, atomCount, name);
    assert.equal((svg.match(/class="lone-pair"/g) || []).length, lonePairs, name);
    assert.equal((svg.match(/class="electron"/g) || []).length, lonePairs * 2, name);
    assert.equal(lines(svg).length, buildLewisGraph(compound).bonds.reduce((sum, bond) => sum + bond.order, 0), name);
    assertBounds(svg);
  }
  const benzene = await resolveCompound('benzene'), svg = renderLewis(benzene);
  assert.equal(labels(svg).length, 12);
  assert.equal(lines(svg).filter(line => line['data-order'] === '2').length, 6);
  assert.equal(lines(svg).length, 15);
});

test('semi-developed formulas keep cyclic graph connectivity and group the exact H total', async () => {
  for (const name of ['benzene','caffeina','5-isobutil decano']) {
    const compound = await resolveCompound(name), svg = renderCondensed(compound);
    const atoms = labels(svg), bondSet = new Set(lines(svg).map(line => line['data-bond']));
    assert.equal(atoms.length, compound.atoms.length, name);
    assert.equal(atoms.reduce((sum, atom) => sum + +atom['data-hydrogens'], 0), compound.atoms.reduce((sum, atom) => sum + atom.hydrogens, 0), name);
    assert.equal(bondSet.size, compound.bonds.length, name);
    assert.ok(svg.includes('Formula semisviluppata'));
    assertBounds(svg);
  }
});

test('unidentified parent chains are never assigned arbitrary carbon numbering', async () => {
  for (const name of ['benzene','caffeina']) {
    const compound = await resolveCompound(name);
    assert.equal((renderTopology(compound, { numbering: true }).match(/<text class="atom-number"/g) || []).length, 0);
  }
});

test('long local chains keep forward bond segments and labels inside view bounds', async () => {
  for (const name of ['decano','eicosano']) {
    const compound = await resolveCompound(name), svg = renderLewis(compound);
    const atomPositions = new Map(labels(svg).map(atom => [+atom['data-atom'], { x: +atom.x, y: +atom.y }]));
    for (const bond of lines(svg)) {
      const [a,b] = bond['data-bond'].split('-').map(Number), p = atomPositions.get(a), q = atomPositions.get(b);
      const dot = (+bond.x2 - +bond.x1) * (q.x - p.x) + (+bond.y2 - +bond.y1) * (q.y - p.y);
      assert.ok(dot > 0, `Reversed or collapsed ${name} bond ${a}-${b}`);
    }
    assertBounds(svg);
  }
});

test('unusual bond orders are dashed placeholders rather than fabricated triple bonds', () => {
  const compound = { name: 'Legame non ordinario', atoms: [{ element:'C', x2d:0, y2d:0, x:0,y:0,z:0, hydrogens:0 }, { element:'C',x2d:1.4,y2d:0,x:1.4,y:0,z:0,hydrogens:0 }], bonds: [{from:0,to:1,order:4}] };
  for (const render of [renderTopology, renderLewis, renderCondensed]) {
    const svg = render(compound), bonds = lines(svg);
    assert.equal(bonds.length, 1);
    assert.ok(bonds[0].class.includes('unknown-bond'));
    assert.equal(bonds[0]['stroke-dasharray'], '4 5');
  }
});

test('names are escaped in SVG markup and rendering leaves chemistry data unchanged', async () => {
  const compound = await resolveCompound('etanolo');
  compound.name = 'Etanolo "<script>alert(1)</script>&';
  const before = JSON.stringify(compound);
  for (const render of [renderTopology, renderLewis, renderCondensed]) {
    const svg = render(compound);
    assert.ok(!svg.includes('<script>'));
    assert.ok(svg.includes('&lt;script&gt;'));
    assert.ok(svg.includes('&quot;'));
    assertBounds(svg);
  }
  assert.equal(JSON.stringify(compound), before);
});

test('3D controls preserve full graph, prefer verified SDF and stop after destruction', async () => {
  const original = globalThis.$3Dmol;
  const calls = [];
  const mock = {
    setBackgroundColor(){}, addModel(...args){ calls.push(['model', ...args]); return {addAtoms(atoms){calls.push(['atoms',atoms]);}}; },
    setStyle(...args){calls.push(['style',...args]);}, render(){}, zoomTo(){}, rotate(){}, getView(){return [0,0,0,-20,0,0,0,1];},
    setView(view){calls.push(['view',view]);}, resize(){calls.push(['resize']);}, spin(...args){calls.push(['spin',...args]);},
    pngURI(){return 'data:image/png;base64,test';}, removeAllModels(){},removeAllLabels(){},removeAllShapes(){}
  };
  globalThis.$3Dmol = {createViewer(){return mock;}};
  const container = {clientWidth:760,clientHeight:380,replaceChildren(){},querySelector(){return null;},querySelectorAll(){return [];}};
  try {
    const compound = await resolveCompound('5-isobutil decano'), viewer = create3DViewer(container,compound);
    assert.equal(viewer.available,true);
    assert.equal(calls.find(call=>call[0]==='model')[2],'sdf');
    viewer.setStyle('spacefill'); viewer.setHydrogens(false); viewer.setSpin(true); viewer.reset(); viewer.resize();
    assert.ok(calls.some(call=>call[0]==='style'&&call[1].elem==='H'&&call[2].sphere.hidden));
    assert.ok(calls.some(call=>call[0]==='spin'&&call[1]==='y'));
    assert.equal(viewer.exportPng(),'data:image/png;base64,test');
    viewer.destroy(); const atDestroy=calls.length; viewer.setSpin(true); viewer.setStyle('stick'); viewer.reset(); viewer.resize();
    assert.equal(calls.length,atDestroy);
    assert.equal(viewer.exportPng(),null);
    calls.length=0;
    const local=await resolveCompound('decano'), estimate=create3DViewer(container,local);
    assert.equal(calls.find(call=>call[0]==='atoms')[1].length,32);
    estimate.destroy();
  } finally {globalThis.$3Dmol=original;}
});
