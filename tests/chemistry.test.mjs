import test from 'node:test';
import assert from 'node:assert/strict';
import { resolveCompound, buildLewisGraph, ChemistryError } from '../dist/chemistry.js';

function neighbors(compound,index) {
  return compound.bonds.flatMap(b=>b.from===index?[b.to]:b.to===index?[b.from]:[]);
}
function verifyValence(compound) {
  const valence=compound.atoms.map(a=>a.hydrogens||0);
  compound.bonds.forEach(b=>{valence[b.from]+=b.order;valence[b.to]+=b.order;});
  compound.atoms.forEach((a,i)=>{if(a.element==='C')assert.equal(valence[i],4,`Carbon ${i} valence`);});
}

test('requested 5-isobutil decano has a decane backbone and the correct isobutyl attachment',async()=>{
  const compound=await resolveCompound('5-isobutil decano');
  assert.equal(compound.formula,'C14H30');
  assert.equal(compound.mass,198.39);
  assert.equal(compound.canonicalName,'5-(2-methylpropyl)decane');
  assert.equal(compound.atoms.length,14);
  assert.equal(compound.bonds.length,13);
  assert.equal(compound.backbone.length,10);
  const parent=compound.backbone[4],main=new Set(compound.backbone);
  const [firstBranch]=neighbors(compound,parent).filter(i=>!main.has(i));
  assert.notEqual(firstBranch,undefined);
  assert.equal(compound.atoms[parent].hydrogens,1);
  assert.equal(compound.atoms[firstBranch].hydrogens,2);
  const [fork]=neighbors(compound,firstBranch).filter(i=>i!==parent);
  assert.equal(compound.atoms[fork].hydrogens,1);
  const methyls=neighbors(compound,fork).filter(i=>i!==firstBranch);
  assert.equal(methyls.length,2);
  methyls.forEach(i=>assert.equal(compound.atoms[i].hydrogens,3));
  assert.equal(compound.condensed,'CH3-(CH2)3-CH(CH2-CH(CH3)2)-(CH2)4-CH3');
  assert.match(compound.conformerSource,/PubChem/);
  assert.match(compound.sdf3d,/ 44 43/);
  assert.ok(compound.atoms.some(atom=>Math.abs(atom.z)>.1));
  verifyValence(compound);
});

test('Italian and English names produce the same branched structure',async()=>{
  for(const name of ['5-isobutil-decano','5-isobutyl decane','5-(2-metilpropil)decano','5-(2-methylpropyl)decane']) {
    const compound=await resolveCompound(name);
    assert.equal(compound.formula,'C14H30');
    assert.match(compound.condensed,/CH2-CH\(CH3\)2/);
    verifyValence(compound);
  }
});

test('multiple locants and mixed substituents preserve carbon valence and hydrogens',async()=>{
  const dimethyl=await resolveCompound('2,3-dimetilesano');
  assert.equal(dimethyl.formula,'C8H18');
  assert.equal(dimethyl.canonicalName,'2,3-dimethylhexane');
  assert.equal(dimethyl.condensed,'CH3-CH(CH3)-CH(CH3)-(CH2)2-CH3');
  verifyValence(dimethyl);
  const mixed=await resolveCompound('3-etil-2-metilpentano');
  assert.equal(mixed.formula,'C8H18');
  assert.equal(mixed.canonicalName,'3-ethyl-2-methylpentane');
  verifyValence(mixed);
  const quaternary=await resolveCompound('2,2-dimetilbutano');
  assert.equal(quaternary.atoms[1].hydrogens,0);
  assert.equal(quaternary.formula,'C6H14');
  verifyValence(quaternary);
  const renumbered=await resolveCompound('5-metilesano');
  assert.equal(renumbered.canonicalName,'2-methylhexane');
  assert.equal(renumbered.backbone[0],5);
  assert.equal(renumbered.condensed,'CH3-CH(CH3)-(CH2)3-CH3');
});

test('iso, sec and tert branches are distinct rooted graphs',async()=>{
  const iso=await resolveCompound('4-isobutilottano');
  const sec=await resolveCompound('4-sec-butilottano');
  const tert=await resolveCompound('4-tert-butilottano');
  assert.equal(iso.formula,'C12H26');
  assert.equal(sec.formula,'C12H26');
  assert.equal(tert.formula,'C12H26');
  assert.equal(iso.atoms[8].hydrogens,2);
  assert.equal(sec.atoms[8].hydrogens,1);
  assert.equal(tert.atoms[8].hydrogens,0);
  [iso,sec,tert].forEach(verifyValence);
});

test('plain chains, rings and unsaturated bonds use the correct hydrogen counts',async()=>{
  const names=[['metano','CH4'],['decano','C10H22'],['eicosano','C20H42'],['cicloesano','C6H12'],['but-2-ene','C4H8'],['2-butene','C4H8'],['pent-1-ino','C5H8'],['3-metiles-1-ene','C7H14']];
  for(const [name,formula]of names) {
    const compound=await resolveCompound(name);
    assert.equal(compound.formula,formula,name);
    verifyValence(compound);
    assert.match(compound.conformerSource,/stimata/);
    assert.ok(compound.notes.some(note=>/stimata/.test(note)));
  }
  assert.equal((await resolveCompound('cicloesano')).condensed,'ciclo-(CH2)6');
});

test('invalid names and impossible locants fail instead of fabricating molecules',async()=>{
  for(const name of ['',null,'<script>benzene</script>','x'.repeat(161)]) {
    await assert.rejects(resolveCompound(name),error=>error instanceof ChemistryError&&error.code==='INVALID_NAME');
  }
  for(const name of ['8-metilpentano','2,2,2-trimetilpropano','2-dimetilesano','2,3-metilesano','but-5-ene']) {
    await assert.rejects(resolveCompound(name),error=>error instanceof ChemistryError&&error.code==='INVALID_NOMENCLATURE');
  }
});

test('explicit Lewis graph conserves atom counts and shows heteroatom lone pairs',async()=>{
  const branched=buildLewisGraph(await resolveCompound('5-isobutil decano'));
  assert.equal(branched.atoms.length,44);
  assert.equal(branched.bonds.length,43);
  assert.equal(branched.atoms.filter(a=>a.element==='H').length,30);
  assert.equal(branched.supported,true);
  const ethanol=buildLewisGraph(await resolveCompound('etanolo'));
  assert.equal(ethanol.atoms.length,9);
  assert.equal(ethanol.atoms.find(a=>a.element==='O').lonePairs,2);
  ethanol.atoms.forEach(a=>['x','y','z','x2d','y2d'].forEach(key=>assert.ok(Number.isFinite(a[key]),key)));
  const caffeine=await resolveCompound('caffeina');
  assert.equal(caffeine.formula,'C8H10N4O2');
  assert.equal(caffeine.condensed,null);
  assert.equal(buildLewisGraph(caffeine).supported,true);
  const benzene=await resolveCompound('benzene');
  assert.equal(benzene.bonds.filter(b=>b.order===2).length,3);
  assert.ok(benzene.notes.some(note=>/delocalizzati/.test(note)));
});

test('cached compounds cannot be corrupted by a caller and aborted lookups stay aborted',async()=>{
  const compound=await resolveCompound('etanolo');
  compound.atoms[0].element='Au';
  assert.notEqual((await resolveCompound('etanolo')).atoms[0].element,'Au');
  const controller=new AbortController();controller.abort();
  await assert.rejects(resolveCompound('etanolo',{signal:controller.signal}),error=>error.name==='AbortError');
});

const mockRecord={
  atoms:{aid:[1,2,3,4,5,6,7,8,9],element:[6,6,8,1,1,1,1,1,1]},
  bonds:{aid1:[1,2,1,1,1,2,2,3],aid2:[2,3,4,5,6,7,8,9],order:[1,1,1,1,1,1,1,1]},
  coords:[{aid:[1,2,3,4,5,6,7,8,9],conformers:[{x:[0,1.3,2.6,-.8,-.8,0,1.3,1.3,3.2],y:[0,.4,0,-.5,.5,1,1.4,-.6,.3]}]}],
};

test('PubChem lookup preserves verified 2D when a 3D conformer is unavailable',async()=>{
  const originalFetch=globalThis.fetch, urls=[];
  globalThis.fetch=async url=>{
    urls.push(url);
    if(url.includes('/property/'))return{ok:true,json:async()=>({PropertyTable:{Properties:[{CID:990001,IUPACName:'ethanol',MolecularFormula:'C2H6O',MolecularWeight:'46.07',ConnectivitySMILES:'CCO'}]}})};
    if(url.includes('record_type=3d'))return{ok:false,status:404};
    return{ok:true,json:async()=>({PC_Compounds:[mockRecord]})};
  };
  try {
    const compound=await resolveCompound('mock-online-ethanol');
    assert.equal(compound.formula,'C2H6O');
    assert.equal(compound.sourceLabel,'PubChem');
    assert.equal(compound.atoms.reduce((sum,a)=>sum+a.hydrogens,0),6);
    assert.equal(compound.sdf3d,undefined);
    assert.match(compound.conformerSource,/stimata/);
    assert.equal(urls.length,3);
    assert.match(compound.sourceUrl,/990001/);
    assert.equal(buildLewisGraph(compound).supported,true);
  } finally {globalThis.fetch=originalFetch;}
});

test('unavailable network and unknown external names have distinct actionable errors',async()=>{
  const originalFetch=globalThis.fetch;
  try {
    globalThis.fetch=async()=>{throw new TypeError('offline');};
    await assert.rejects(resolveCompound('unknown-online-structure'),error=>error.code==='NETWORK');
    globalThis.fetch=async()=>({ok:false,status:404});
    await assert.rejects(resolveCompound('not-a-chemical-name'),error=>error.code==='NOT_FOUND');
    // Local structures remain available after failed external requests.
    assert.equal((await resolveCompound('propano')).formula,'C3H8');
  } finally {globalThis.fetch=originalFetch;}
});

test('ambiguous names do not silently select one of several structures',async()=>{
  const originalFetch=globalThis.fetch;
  globalThis.fetch=async()=>({ok:true,json:async()=>({PropertyTable:{Properties:[{CID:1},{CID:2}]}})});
  try {await assert.rejects(resolveCompound('ambiguous-test-mixture'),error=>error.code==='AMBIGUOUS_NAME');}
  finally {globalThis.fetch=originalFetch;}
});

test('external carbonyl branches keep bond order and charged groups keep formal charge',async()=>{
  const originalFetch=globalThis.fetch;
  const ketone={
    atoms:{aid:[1,2,3,4,5,6,7,8,9,10,11,12,13],element:[8,6,6,6,6,1,1,1,1,1,1,1,1]},
    bonds:{aid1:[1,2,2,3,3,3,4,4,4,5,5,5],aid2:[2,3,4,5,6,7,8,9,10,11,12,13],order:[2,1,1,1,1,1,1,1,1,1,1,1]},
    coords:[{aid:[1,2,3,4,5,6,7,8,9,10,11,12,13],conformers:[{x:[0,0,1.3,-1.3,2.6,1,1,-1,-1,-1,3,3,3],y:[1.4,0,-.4,-.4,0,1,-1,1,-1,0,1,0,-1]}]}],
  };
  const ammonium={
    atoms:{aid:[1,2,3,4,5],element:[7,1,1,1,1],charge:[{aid:1,value:1}]},
    bonds:{aid1:[1,1,1,1],aid2:[2,3,4,5],order:[1,1,1,1]},
    coords:[{aid:[1,2,3,4,5],conformers:[{x:[0,1,-1,0,0],y:[0,0,0,1,-1]}]}],
  };
  let current='ketone';
  globalThis.fetch=async url=>{
    if(url.includes('/property/'))return{ok:true,json:async()=>({PropertyTable:{Properties:[current==='ketone'?{CID:990002,IUPACName:'butan-2-one',MolecularFormula:'C4H8O',MolecularWeight:'72.11'}:{CID:990003,IUPACName:'ammonium',MolecularFormula:'H4N+',MolecularWeight:'18.04'}]}})};
    if(url.includes('record_type=3d'))return{ok:false,status:404};
    return{ok:true,json:async()=>({PC_Compounds:[current==='ketone'?ketone:ammonium]})};
  };
  try {
    const compound=await resolveCompound('mock-carbonyl-branch');
    assert.match(compound.condensed,/C\(=O\)/);
    assert.equal(compound.backbone,undefined,'heteroatom chain is not exposed as numbered carbon parent');
    current='ammonium';
    const charged=await resolveCompound('mock-charged-ammonium');
    assert.equal(charged.condensed,'[NH4+]');
    const lewis=buildLewisGraph(charged);
    assert.equal(lewis.atoms[0].charge,1);
    assert.equal(lewis.atoms[0].lonePairs,0);
    assert.equal(lewis.supported,true);
  } finally {globalThis.fetch=originalFetch;}
});

test('cyclic unsaturated nomenclature outside the local subset reaches verified lookup',async()=>{
  const originalFetch=globalThis.fetch;let called=false;
  globalThis.fetch=async()=>{called=true;return{ok:false,status:404};};
  try {
    await assert.rejects(resolveCompound('cyclohexene'),error=>error.code==='NOT_FOUND');
    assert.equal(called,true);
  } finally {globalThis.fetch=originalFetch;}
});
