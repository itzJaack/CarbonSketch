/**
 * CarbonSketch chemistry: local nomenclature and public PubChem records.
 * PubChem PUG REST: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
 * Bundled reference coordinates were retrieved from PubChem on 2026-10-07.
 * Locally generated conformers are always explicitly labelled as estimates.
 */

const API = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
const ATOMIC_WEIGHTS = { H:1.008, B:10.81, C:12.011, N:14.007, O:15.999, F:18.998, Si:28.085, P:30.974, S:32.06, Cl:35.45, Br:79.904, I:126.904 };
const ELEMENTS = ('? H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split(' ');
const VALENCE_ELECTRONS = { H:1, B:3, C:4, N:5, O:6, F:7, Si:4, P:5, S:6, Cl:7, Br:7, I:7 };
const cache = new Map();

export class ChemistryError extends Error {
  constructor(message, code='UNKNOWN') { super(message); this.name='ChemistryError'; this.code=code; }
}

function cleanName(name) {
  if (typeof name !== 'string') throw new ChemistryError('Inserisci il nome di un composto.', 'INVALID_NAME');
  const value = name.normalize('NFKC').trim();
  if (!value || value.length > 160 || /[<>\u0000-\u001f]/.test(value)) throw new ChemistryError('Inserisci un nome chimico valido, fino a 160 caratteri.', 'INVALID_NAME');
  return value;
}
function normalize(name) {
  return name.toLowerCase().normalize('NFD').replace(/[\u0300-\u036f]/g,'').replace(/[‐‑–—−]/g,'-').replace(/\s+/g,'');
}
function abortCheck(signal) { if (signal?.aborted) throw new DOMException('Ricerca annullata.', 'AbortError'); }
function clone(value) { return structuredClone(value); }
function adjacency(compound) {
  const graph = compound.atoms.map(()=>[]);
  compound.bonds.forEach(b => { graph[b.from].push({to:b.to,order:b.order}); graph[b.to].push({to:b.from,order:b.order}); });
  return graph;
}
function formulaFromGraph(atoms) {
  const counts = {};
  atoms.forEach(a => { counts[a.element]=(counts[a.element]||0)+1; counts.H=(counts.H||0)+(a.hydrogens||0); });
  const order = counts.C ? ['C','H',...Object.keys(counts).filter(e=>e!=='C'&&e!=='H').sort()] : Object.keys(counts).sort();
  return order.filter(e=>counts[e]).map(e=>e+(counts[e]===1?'':counts[e])).join('');
}
function graphMass(atoms) { return atoms.reduce((sum,a)=>sum+(ATOMIC_WEIGHTS[a.element]||0)+(a.hydrogens||0)*ATOMIC_WEIGHTS.H,0); }

const CHAIN_ROOTS = [
  [20,['eicos','icos']], [19,['nonadec']], [18,['octadec','ottadec']], [17,['heptadec','eptadec']], [16,['hexadec','esadec']],
  [15,['pentadec']], [14,['tetradec']], [13,['tridec']], [12,['dodec']], [11,['undec']], [10,['dec']],
  [9,['non']], [8,['oct','ott']], [7,['hept','ept']], [6,['hex','es']], [5,['pent']], [4,['but']], [3,['prop']], [2,['eth','et']], [1,['meth','met']],
];
const ENGLISH_ROOTS = {1:'meth',2:'eth',3:'prop',4:'but',5:'pent',6:'hex',7:'hept',8:'oct',9:'non',10:'dec',11:'undec',12:'dodec',13:'tridec',14:'tetradec',15:'pentadec',16:'hexadec',17:'heptadec',18:'octadec',19:'nonadec',20:'eicos'};
// Match the longest textual root first: metano/methane must not match etano/ethane.
const SORTED_ROOTS=CHAIN_ROOTS.flatMap(([length,roots])=>roots.map(root=>({length,root}))).sort((a,b)=>b.root.length-a.root.length);
// Edges are rooted at atom 0: isobutyl means -CH2-CH(CH3)2, never tert-butyl.
const BRANCHES = [
  {names:['(2-methylpropyl)','(2-metilpropil)','isobutyl','isobutil'],canonical:'(2-methylpropyl)',count:4,edges:[[0,1],[1,2],[1,3]]},
  {names:['(1-methylethyl)','(1-metiletil)','isopropyl','isopropil'],canonical:'(propan-2-yl)',count:3,edges:[[0,1],[0,2]]},
  {names:['tert-butyl','tert-butil','t-butyl','t-butil','terz-butil','(1,1-dimethylethyl)','(1,1-dimetiletil)'],canonical:'tert-butyl',count:4,edges:[[0,1],[0,2],[0,3]]},
  {names:['sec-butyl','sec-butil','s-butyl','s-butil','(1-methylpropyl)','(1-metilpropil)'],canonical:'(butan-2-yl)',count:4,edges:[[0,1],[0,2],[2,3]]},
  ...[[1,'methyl','metil'],[2,'ethyl','etil'],[3,'propyl','propil'],[4,'butyl','butil'],[5,'pentyl','pentil'],[6,'hexyl','esil']].map(([count,en,it])=>({names:[`n-${en}`,`n-${it}`,en,it],canonical:en,count,edges:Array.from({length:count-1},(_,i)=>[i,i+1])})),
];
const BRANCH_NAMES = BRANCHES.flatMap(branch=>branch.names.map(token=>({token,branch}))).sort((a,b)=>b.token.length-a.token.length);
const MULTIPLIERS = [{token:'deca',count:10},{token:'nona',count:9},{token:'octa',count:8},{token:'otta',count:8},{token:'hepta',count:7},{token:'epta',count:7},{token:'tetra',count:4},{token:'penta',count:5},{token:'hexa',count:6},{token:'esa',count:6},{token:'tri',count:3},{token:'di',count:2}];
const ITALIAN_MULTIPLIERS={2:'di',3:'tri',4:'tetra',5:'penta',6:'esa',7:'epta',8:'otta',9:'nona',10:'deca'};
const ENGLISH_MULTIPLIERS={...ITALIAN_MULTIPLIERS,6:'hexa',7:'hepta',8:'octa'};
const ITALIAN_ROOTS={1:'met',2:'et',3:'prop',4:'but',5:'pent',6:'es',7:'ept',8:'ott',9:'non',10:'dec',11:'undec',12:'dodec',13:'tridec',14:'tetradec',15:'pentadec',16:'esadec',17:'eptadec',18:'ottadec',19:'nonadec',20:'eicos'};
const ITALIAN_BRANCHES={'methyl':'metil','ethyl':'etil','propyl':'propil','butyl':'butil','pentyl':'pentil','hexyl':'esil','(2-methylpropyl)':'isobutil','(propan-2-yl)':'isopropil','(butan-2-yl)':'sec-butil','tert-butyl':'tert-butil'};
const alphabeticalBranch=label=>label.replace(/^tert-/,'').replace(/[()0-9,-]/g,'');

function findParent(normalized) {
  for (const {length,root} of SORTED_ROOTS) {
    for (const suffix of ['ane','ano','ene','eno','yne','ino']) {
      const endings = [root+suffix];
      for (const ending of endings) {
        if (!normalized.endsWith(ending)) continue;
        let prefix=normalized.slice(0,-ending.length), cycle=false, position=1;
        if (prefix.endsWith('cyclo')) { cycle=true; prefix=prefix.slice(0,-5); }
        if (prefix.endsWith('ciclo')) { cycle=true; prefix=prefix.slice(0,-5); }
        const order = suffix==='ene'||suffix==='eno'?2:suffix==='yne'||suffix==='ino'?3:1;
        if (order>1) {
          const match=prefix.match(/(?:^|-)(\d+)-?$/);
          if (match) { position=Number(match[1]); prefix=prefix.slice(0,match.index); }
        }
        return {length,prefix:prefix.replace(/-$/,''),cycle,order,position};
      }
      // Preferred spelling: hex-2-ene / es-2-ene / pent-1-ino.
      const match=normalized.match(new RegExp(`^(.*)${root}-(\\d+)-(ene|eno|yne|ino)$`));
      if (match) return {length,prefix:match[1].replace(/-$/,''),cycle:false,order:/ene|eno/.test(match[3])?2:3,position:Number(match[2])};
    }
  }
  return null;
}
function parseBranches(prefix,parentLength) {
  const result=[];
  let rest=prefix;
  while(rest) {
    // Locants may be comma-separated, dash-separated, and adjacent to the group.
    // This is restricted to digit lists followed by a recognized branch token.
    const locants=rest.match(/^(\d+(?:[,-]\d+)*)(?:-)?/);
    if(!locants) return null;
    const positions=locants[1].split(/[,-]/).map(Number); rest=rest.slice(locants[0].length);
    let multiplier=1;
    const multi=MULTIPLIERS.find(m=>rest.startsWith(m.token)&&BRANCH_NAMES.some(b=>rest.slice(m.token.length).replace(/^-/,'').startsWith(b.token)));
    if(multi) {multiplier=multi.count; rest=rest.slice(multi.token.length).replace(/^-/,'');}
    const match=BRANCH_NAMES.find(b=>rest.startsWith(b.token));
    if(!match) return null;
    rest=rest.slice(match.token.length).replace(/^-/,'');
    if(multiplier!==positions.length) throw new ChemistryError(`Il gruppo ${multi?.token||match.token} richiede ${multiplier} locanti, ma ne hai indicati ${positions.length}.`,'MULTIPLIER_MISMATCH');
    positions.forEach(position=>{
      if(position<1||position>parentLength) throw new ChemistryError(`Il locante ${position} è fuori dalla catena di ${parentLength} carboni.`,'INVALID_LOCANT');
      result.push({position,...match.branch});
    });
  }
  return result;
}

function chainFragment(compound,backbone) {
  const graph=adjacency(compound),mainSet=new Set(backbone),seen=new Set(backbone);
  const atomText=index=>{const a=compound.atoms[index];const symbol=a.element+(a.hydrogens?'H'+(a.hydrogens===1?'':a.hydrogens):'');return a.charge?'['+symbol+(Math.abs(a.charge)===1?'':Math.abs(a.charge))+(a.charge>0?'+':'-')+']':symbol;};
  function branchText(index,parent) {
    seen.add(index);
    const children=graph[index].filter(n=>n.to!==parent&&!seen.has(n.to));
    const head=atomText(index);
    // Identical terminal groups are compacted, e.g. CH(CH3)2.
    if(children.length>1&&children.every(n=>graph[n.to].length===1&&atomText(n.to)===atomText(children[0].to))) {
      children.forEach(n=>seen.add(n.to));return head+'('+atomText(children[0].to)+')'+children.length;
    }
    return head+children.map((n,i)=>{
      const edge=n.order===2?'=':n.order===3?'≡':'';
      const text=branchText(n.to,index);
      return children.length===1?edge+text:'('+edge+text+')';
    }).join('');
  }
  const segments=backbone.map((index,i)=>{
    const branches=graph[index].filter(n=>!mainSet.has(n.to));
    let text=atomText(index)+branches.map(n=>'('+(n.order===2?'=':n.order===3?'≡':'')+branchText(n.to,index)+')').join('');
    const edge=i<backbone.length-1?graph[index].find(n=>n.to===backbone[i+1])?.order:0;
    return {text,edge};
  });
  // Repeated unbranched methylene groups become (CH2)n in condensed form.
  const output=[];
  for(let i=0;i<segments.length;i++) {
    let end=i;
    while(segments[i].text==='CH2'&&segments[i].edge===1&&end+1<segments.length&&segments[end+1].text==='CH2'&&segments[end+1].edge===1) end++;
    const count=end-i+1;
    output.push(count>1?`(CH2)${count}`:segments[i].text);
    const edge=segments[end].edge;
    if(edge) output.push(edge===2?'=':edge===3?'≡':'-');
    i=end;
  }
  return output.join('');
}

function parseLocal(name) {
  const parent=findParent(normalize(name));
  if(!parent) return null;
  const branches=parseBranches(parent.prefix,parent.length);
  if(!branches) return null;
  const {length,cycle,order,position}=parent;
  if(cycle&&length<3) throw new ChemistryError('Una catena ciclica deve contenere almeno tre carboni.','INVALID_NOMENCLATURE');
  if(cycle&&order>1)return null; // Valid nomenclature outside the local subset goes to PubChem.
  if(order>1&&(length<2||position<1||position>=length)) throw new ChemistryError('Il locante del legame multiplo è fuori dalla catena.','INVALID_LOCANT');
  const atoms=[],bonds=[],backbone=Array.from({length},(_,i)=>i);
  for(let i=0;i<length;i++) {
    const angle=2*Math.PI*i/length;
    const x2d=cycle?Math.cos(angle)*length*.24:(i-(length-1)/2)*1.3;
    const y2d=cycle?Math.sin(angle)*length*.24:(i%2?.39:-.39);
    atoms.push({element:'C',x:x2d,y:y2d,z:cycle&&length===6?(i%2?.3:-.3):0,x2d,y2d,charge:0,hydrogens:0});
    if(i) bonds.push({from:i-1,to:i,order:i===position?order:1});
  }
  if(cycle) bonds.push({from:length-1,to:0,order:1});
  const attachmentCount={};
  branches.forEach(branch=>{
    const base=branch.position-1, offset=atoms.length;
    const slot=attachmentCount[base]||0; attachmentCount[base]=slot+1;
    const direction=slot%2===0?1:-1;
    const branchPositions=[{x2d:atoms[base].x2d,y2d:atoms[base].y2d+direction*1.6,x:atoms[base].x+.16,y:atoms[base].y+direction*.64,z:direction*1.34}];
    branch.edges.forEach(([a,b])=>{
      const p=branchPositions[a];
      const siblings=branch.edges.filter(edge=>edge[0]===a);
      const side=siblings.length>1?(siblings.findIndex(edge=>edge[1]===b)%2===0?-1:1):(b%2===0?-.5:.5);
      branchPositions[b]={x2d:p.x2d+side*1.15,y2d:p.y2d+direction*1.28,x:p.x+side*1.2,y:p.y+direction*1.1,z:p.z+(b%2?-.25:.25)};
    });
    branchPositions.forEach(p=>atoms.push({element:'C',...p,charge:0,hydrogens:0}));
    bonds.push({from:base,to:offset,order:1});
    branch.edges.forEach(([a,b])=>bonds.push({from:offset+a,to:offset+b,order:1}));
  });
  const degree=atoms.map(()=>0); bonds.forEach(b=>{degree[b.from]+=b.order;degree[b.to]+=b.order;});
  atoms.forEach((a,i)=>{
    if(degree[i]>4) throw new ChemistryError(`Il carbonio ${i+1} avrebbe ${degree[i]} legami: un carbonio neutro ne può formare quattro. Controlla i sostituenti e i locanti.`,'INVALID_VALENCE');
    a.hydrogens=4-degree[i];
  });
  const family=cycle?'Cicloalcano':order===1?(branches.length?'Alcano ramificato':'Alcano'):order===2?'Alchene':'Alchino';
  // Normalize a simple acyclic parent direction using the lowest locant set.
  const compareLocants=(a,b)=>{for(let i=0;i<Math.max(a.length,b.length);i++){if((a[i]??0)!==(b[i]??0))return(a[i]??0)-(b[i]??0);}return 0;};
  const forwardLocants=branches.map(b=>b.position).sort((a,b)=>a-b);
  const reverseLocants=branches.map(b=>length+1-b.position).sort((a,b)=>a-b);
  const reverse=!cycle&&(order>1?(length-position<position||length-position===position&&compareLocants(reverseLocants,forwardLocants)<0):compareLocants(reverseLocants,forwardLocants)<0);
  if(reverse)backbone.reverse();
  const grouped=new Map();branches.forEach(b=>{const key=b.canonical;grouped.set(key,[...(grouped.get(key)||[]),reverse?length+1-b.position:b.position]);});
  const groups=[...grouped].sort((a,b)=>alphabeticalBranch(a[0]).localeCompare(alphabeticalBranch(b[0])));
  const canonicalPrefix=groups.map(([label,locants])=>`${locants.sort((a,b)=>a-b).join(',')}-${locants.length===1?'':(ENGLISH_MULTIPLIERS[locants.length]||'')}${label}`).join('-');
  const canonicalParent=(cycle?'cyclo':'')+ENGLISH_ROOTS[length]+(order===1?'ane':`-${reverse?length-position:position}-${order===2?'ene':'yne'}`);
  const italianPrefix=groups.map(([label,locants])=>`${locants.join(',')}-${locants.length===1?'':(ITALIAN_MULTIPLIERS[locants.length]||'')}${ITALIAN_BRANCHES[label]||label}`).join('-');
  const normalizedName=italianPrefix+(cycle?'ciclo':'')+ITALIAN_ROOTS[length]+(order===1?'ano':`-${reverse?length-position:position}-${order===2?'ene':'ino'}`);
  const compound={name,canonicalName:(canonicalPrefix||'')+canonicalParent,canonicalNameVerified:false,namingStatus:'interpreted',normalizedName,inputInterpretation:`Interpretato come ${normalizedName}.`,formula:formulaFromGraph(atoms),mass:Number(graphMass(atoms).toFixed(3)),family,sourceLabel:'Nomenclatura locale',sourceUrl:'https://iupac.org/what-we-do/nomenclature/',atoms,bonds,backbone,conformerSource:'Geometria 3D stimata · non ottimizzata',notes:['La geometria 3D è una rappresentazione didattica stimata, non un conformero calcolato con un metodo energetico.']};
  compound.condensed=cycle?(branches.length?`ciclo-(${chainFragment(compound,backbone)})`:`ciclo-(CH2)${length}`):chainFragment(compound,backbone);
  if(!cycle&&longestPath(compound)?.length>length)compound.notes.push('Nome non preferito: la struttura è interpretabile, ma la catena principale indicata non è la più lunga.');
  if(cycle&&branches.length)compound.notes.push('Per gli anelli sostituiti il nome locale è una traduzione del nome inserito; la numerazione IUPAC dell’anello non viene ricostruita.');
  return compound;
}

function readPubChem(record) {
  if(!record?.atoms?.aid?.length) throw new ChemistryError('PubChem non ha restituito una struttura utilizzabile.','INVALID_STRUCTURE');
  const ids=record.atoms.aid,charges=new Map((record.atoms.charge||[]).map(item=>[item.aid,item.value]));
  const coords=record.coords?.[0],conformer=coords?.conformers?.[0];
  const coordMap=new Map((coords?.aid||[]).map((id,i)=>[id,{x:conformer?.x?.[i]??0,y:conformer?.y?.[i]??0,z:conformer?.z?.[i]??0}]));
  const fullAtoms=ids.map((id,i)=>({aid:id,element:ELEMENTS[record.atoms.element[i]]||'?',charge:charges.get(id)||0,...coordMap.get(id)}));
  const fullMap=new Map(fullAtoms.map((atom,i)=>[atom.aid,i]));
  const fullBonds=(record.bonds?.aid1||[]).map((id,i)=>({from:fullMap.get(id),to:fullMap.get(record.bonds.aid2[i]),order:record.bonds.order[i]}));
  const heavy=fullAtoms.filter(a=>a.element!=='H');
  if(!heavy.length) throw new ChemistryError('Le strutture composte soltanto da idrogeno non sono supportate dal visualizzatore.','UNSUPPORTED_STRUCTURE');
  if(heavy.length>100) throw new ChemistryError('La struttura supera il limite di 100 atomi non idrogeno del visualizzatore.','STRUCTURE_TOO_LARGE');
  const heavyMap=new Map(heavy.map((a,i)=>[a.aid,i]));
  const atoms=heavy.map(a=>({...a,x2d:a.x,y2d:a.y,hydrogens:0}));
  const bonds=[];
  fullBonds.forEach(b=>{
    const a=fullAtoms[b.from],c=fullAtoms[b.to];
    if(a.element==='H'&&c.element!=='H') atoms[heavyMap.get(c.aid)].hydrogens++;
    else if(c.element==='H'&&a.element!=='H') atoms[heavyMap.get(a.aid)].hydrogens++;
    else if(a.element!=='H'&&c.element!=='H') bonds.push({from:heavyMap.get(a.aid),to:heavyMap.get(c.aid),order:b.order});
  });
  return {atoms,bonds,fullAtoms,fullBonds,is3d:Boolean(conformer?.z?.length)};
}
function sdfFromRecord(data,title='PubChem 3D conformer') {
  const {fullAtoms:atoms,fullBonds:bonds}=data;
  const pad=(number,width)=>String(number).padStart(width,' ');
  const lines=[title,'  CarbonSketch        3D','',`${pad(atoms.length,3)}${pad(bonds.length,3)}  0  0  0  0            999 V2000`];
  atoms.forEach(a=>lines.push(`${a.x.toFixed(4).padStart(10)}${a.y.toFixed(4).padStart(10)}${a.z.toFixed(4).padStart(10)} ${a.element.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`));
  bonds.forEach(b=>lines.push(`${pad(b.from+1,3)}${pad(b.to+1,3)}${pad(b.order,3)}  0  0  0  0`));
  const charged=atoms.map((a,i)=>({index:i+1,charge:a.charge})).filter(a=>a.charge);
  for(let i=0;i<charged.length;i+=8) { const chunk=charged.slice(i,i+8);lines.push(`M  CHG${pad(chunk.length,3)}`+chunk.map(a=>`${pad(a.index,4)}${pad(a.charge,4)}`).join('')); }
  lines.push('M  END','$$$$');
  return lines.join('\n');
}
function merge3d(compound,data) {
  if(!data?.is3d) return false;
  const byId=new Map(data.atoms.map(a=>[a.aid,a]));
  if(compound.atoms.some(a=>byId.get(a.aid)?.element!==a.element)) return false;
  compound.atoms.forEach(a=>{const other=byId.get(a.aid);a.x=other.x;a.y=other.y;a.z=other.z;});
  compound.sdf3d=sdfFromRecord(data,compound.canonicalName);
  compound.conformerSource='Conformero 3D PubChem';
  return true;
}
function longestPath(compound) {
  const graph=adjacency(compound);if(compound.bonds.length!==compound.atoms.length-1) return null;
  function farthest(start) {
    const stack=[{node:start,path:[start]}];let best=[start];
    while(stack.length) {const {node,path}=stack.pop();if(path.length>best.length)best=path;graph[node].forEach(n=>{if(!path.includes(n.to))stack.push({node:n.to,path:[...path,n.to]});});}return best;
  }
  const endpoint=farthest(0).at(-1);return farthest(endpoint);
}
function condensedForExternal(compound) {
  const backbone=longestPath(compound);
  if(backbone&&compound.atoms.every(a=>['C','N','O','F','Cl','Br','I'].includes(a.element))) {
    if(compound.atoms.every(a=>a.element==='C'))compound.backbone=backbone;
    return chainFragment(compound,backbone);
  }
  compound.notes.push('Per questa struttura ciclica o complessa è mostrata una formula semisviluppata con tutti i collegamenti, al posto di una condensata lineare.');
  return null;
}
function inferFamily(atoms,bonds) {
  if(atoms.some(a=>!['C','H'].includes(a.element))) return 'Composto organico';
  if(bonds.some(b=>b.order===3))return 'Idrocarburo insaturo';
  if(bonds.some(b=>b.order===2))return 'Idrocarburo insaturo';
  return bonds.length===atoms.length?'Cicloalcano':'Alcano';
}

let requestQueue=Promise.resolve(), lastRequest=0;
async function request(path,signal) {
  abortCheck(signal);
  const turn=requestQueue.then(async()=>{
    const wait=Math.max(0,250-(Date.now()-lastRequest));if(wait)await new Promise(resolve=>setTimeout(resolve,wait));
    abortCheck(signal);lastRequest=Date.now();
  });
  requestQueue=turn.catch(()=>{});await turn;
  const controller=new AbortController(),onAbort=()=>controller.abort();
  signal?.addEventListener('abort',onAbort,{once:true});
  const timer=setTimeout(()=>controller.abort(),12000);
  try {
    const response=await fetch(API+path,{signal:controller.signal,headers:{Accept:'application/json'}});
    if(!response.ok)throw new ChemistryError(response.status===404?'Composto non trovato in PubChem. Prova il nome IUPAC in inglese.':response.status===503?'PubChem è temporaneamente occupato. Riprova tra poco.':'PubChem non è disponibile in questo momento.',response.status===404?'NOT_FOUND':'NETWORK');
    return await response.json();
  } catch(error) {
    abortCheck(signal);
    if(error instanceof ChemistryError)throw error;
    throw new ChemistryError('Impossibile contattare PubChem. Controlla la connessione oppure prova uno degli esempi disponibili offline.','NETWORK');
  } finally {clearTimeout(timer);signal?.removeEventListener('abort',onAbort);}
}

const ITALIAN_ALIASES={'acidoacetilsalicilico':'aspirin','aspirina':'aspirin','glucosio':'glucose','fruttosio':'fructose','saccarosio':'sucrose','paracetamolo':'acetaminophen','ibuprofene':'ibuprofen','ammoniaca':'ammonia','anidridecarbonica':'carbon dioxide','clorurodisodio':'sodium chloride','toluene':'toluene','fenolo':'phenol','anilina':'aniline','formaldeide':'formaldehyde','etilene':'ethylene','acetilene':'acetylene','glicerolo':'glycerol','glicerina':'glycerol','acidoformico':'formic acid','acidolattico':'lactic acid','acidocitrico':'citric acid'};

async function lookupPubChem(name,signal) {
  const query=ITALIAN_ALIASES[normalize(name)]||name;
  const result=await request(`/compound/name/${encodeURIComponent(query)}/property/MolecularFormula,MolecularWeight,IUPACName,ConnectivitySMILES/JSON`,signal);
  const properties=result.PropertyTable?.Properties;
  if(!properties?.length)throw new ChemistryError('Composto non trovato. Prova un nome IUPAC più preciso.','NOT_FOUND');
  if(properties.length>1)throw new ChemistryError('Il nome identifica più composti. Specifica il nome IUPAC completo.','AMBIGUOUS_NAME');
  const p=properties[0];
  const [two,three]=await Promise.allSettled([
    request(`/compound/cid/${p.CID}/JSON?record_type=2d`,signal),
    request(`/compound/cid/${p.CID}/JSON?record_type=3d`,signal),
  ]);
  abortCheck(signal);
  if(two.status!=='fulfilled')throw two.reason;
  const graph=readPubChem(two.value.PC_Compounds?.[0]);
  const compound={name,canonicalName:p.IUPACName||query,canonicalNameVerified:Boolean(p.IUPACName),namingStatus:'verified',normalizedName:p.IUPACName||query,formula:p.MolecularFormula,mass:Number(p.MolecularWeight),family:inferFamily(graph.atoms,graph.bonds),sourceLabel:'PubChem',sourceUrl:`https://pubchem.ncbi.nlm.nih.gov/compound/${p.CID}`,smiles:p.ConnectivitySMILES||p.SMILES,atoms:graph.atoms,bonds:graph.bonds,conformerSource:'Geometria 3D stimata · conformero non disponibile',notes:[]};
  let merged=false;
  if(three.status==='fulfilled') {try{merged=merge3d(compound,readPubChem(three.value.PC_Compounds?.[0]));}catch{/* The verified 2D structure remains usable when 3D is unavailable. */}}
  if(!merged) {
    // A plane is an honest preview when no real conformer is available.
    compound.atoms.forEach(a=>a.z=0);
    compound.notes.push('PubChem non ha fornito un conformero 3D: la vista mostra una geometria didattica planare stimata.');
  }
  if(graph.atoms.some(a=>a.element!=='C'&&a.element!=='H')&&graph.atoms.every(a=>a.element!=='C'))compound.family='Composto inorganico';
  compound.condensed=condensedForExternal(compound);
  if(compound.bonds.some(b=>![1,2,3].includes(b.order)))compound.notes.push('Questa struttura contiene legami non ordinari; la rappresentazione di Lewis completa non è disponibile.');
  return compound;
}

/** Resolve an Italian/English alkane name locally, or a verified PubChem name online. */
export async function resolveCompound(input,{signal}={}) {
  const name=cleanName(input);abortCheck(signal);
  const key=normalize(name);
  if(cache.has(key)) {const result=clone(cache.get(key));result.name=name;return result;}
  const bundled=resolveBundled(name,key);
  const compound=bundled||parseLocal(name)||await lookupPubChem(name,signal);
  abortCheck(signal);cache.set(key,clone(compound));
  return compound;
}

function graphValences(element,charge) {
  if(!Number.isInteger(charge)||Math.abs(charge)>1)return null;
  if(element==='H')return charge===0?[1]:[0];
  if(element==='C'||element==='Si')return charge===0?[4]:[3];
  if(element==='B')return charge===0?[3]:charge===-1?[4]:[2];
  if(element==='N')return charge===0?[3]:charge===1?[4]:[2];
  if(element==='O')return charge===0?[2]:charge===1?[3]:[1];
  if(['F','Cl','Br','I'].includes(element))return charge===0?[1]:charge===-1?[0]:[2];
  if(element==='P')return charge===0?[3,5]:charge===1?[4,6]:[2,4,6];
  if(element==='S')return charge===0?[2,4,6]:charge===1?[3,5]:[1,3,5];
  return null;
}

function rootedSignature(graph,node,parent=-1) {
  return '('+graph[node].filter(n=>n.to!==parent).map(n=>rootedSignature(graph,n.to,node)).sort().join('')+')';
}
function referenceSignature(branch) {
  const graph=Array.from({length:branch.count},()=>[]);
  branch.edges.forEach(([a,b])=>{graph[a].push({to:b});graph[b].push({to:a});});
  return rootedSignature(graph,0);
}
const BRANCH_SIGNATURES=BRANCHES.map(branch=>({signature:referenceSignature(branch),branch}));
function lexicographicNumbers(a,b) {
  for(let i=0;i<Math.max(a.length,b.length);i++){const difference=(a[i]??0)-(b[i]??0);if(difference)return difference;}
  return 0;
}
function interpretedTreeName(compound) {
  if(!compound.atoms.every(a=>a.element==='C'&&a.charge===0)||compound.bonds.some(b=>b.order!==1)||compound.bonds.length!==compound.atoms.length-1)return null;
  const graph=adjacency(compound),leaves=graph.flatMap((neighbors,i)=>neighbors.length<=1?[i]:[]),candidates=[];
  let maximumLength=0;
  for(const start of leaves) {
    const stack=[{node:start,parent:-1,path:[start]}];
    while(stack.length) {
      const {node,parent,path}=stack.pop();
      if(graph[node].length<=1) {
        if(path.length>maximumLength){maximumLength=path.length;candidates.length=0;}
        if(path.length===maximumLength)candidates.push(path);
      }
      graph[node].forEach(n=>{if(n.to!==parent)stack.push({node:n.to,parent:node,path:[...path,n.to]});});
    }
  }
  if(!ENGLISH_ROOTS[maximumLength])return null;
  const interpreted=[];
  for(const path of candidates) {
    const mainSet=new Set(path),groups=new Map();let valid=true;
    path.forEach((node,index)=>graph[node].filter(n=>!mainSet.has(n.to)).forEach(n=>{
      const signature=rootedSignature(graph,n.to,node);
      const recognized=BRANCH_SIGNATURES.find(item=>item.signature===signature)?.branch;
      if(!recognized){valid=false;return;}
      const list=groups.get(recognized.canonical)||[];list.push(index+1);groups.set(recognized.canonical,list);
    }));
    if(!valid||[...groups.values()].some(locants=>locants.length>1&&!ENGLISH_MULTIPLIERS[locants.length]))continue;
    const sorted=[...groups].sort((a,b)=>alphabeticalBranch(a[0]).localeCompare(alphabeticalBranch(b[0])));
    const allLocants=sorted.flatMap(([,locants])=>locants).sort((a,b)=>a-b);
    const enPrefix=sorted.map(([label,locants])=>locants.join(',')+'-'+(locants.length===1?'':ENGLISH_MULTIPLIERS[locants.length])+label).join('-');
    const itPrefix=sorted.map(([label,locants])=>locants.join(',')+'-'+(locants.length===1?'':ITALIAN_MULTIPLIERS[locants.length])+(ITALIAN_BRANCHES[label]||label)).join('-');
    interpreted.push({backbone:path,canonicalName:enPrefix+ENGLISH_ROOTS[path.length]+'ane',normalizedName:itPrefix+ITALIAN_ROOTS[path.length]+'ano',branchCount:allLocants.length,locants:allLocants,alphabeticalLocants:sorted.flatMap(([,locants])=>locants)});
  }
  interpreted.sort((a,b)=>b.branchCount-a.branchCount||lexicographicNumbers(a.locants,b.locants)||lexicographicNumbers(a.alphabeticalLocants,b.alphabeticalLocants));
  return interpreted[0]||null;
}

/**
 * Derive a compound from an explicitly drawn graph, without mutating it.
 * H counts are exact when supplied; otherwise the lowest supported valence is used.
 * Only simple connected graphs and the supported closed-shell element valences
 * are accepted. Metadata and atom groups are retained; imported SDF provenance
 * is discarded because drawn coordinates are always estimated.
 */
export function compoundFromGraph(input,options={}) {
  if(!input||!Array.isArray(input.atoms)||!Array.isArray(input.bonds)||!input.atoms.length)throw new ChemistryError('Disegna almeno un atomo e fornisci una lista di legami.','GRAPH_EMPTY');
  if(input.atoms.length>100)throw new ChemistryError('Il disegno può contenere al massimo 100 atomi.','STRUCTURE_TOO_LARGE');
  const graphInput=clone(input),atoms=graphInput.atoms,bonds=graphInput.bonds;
  const bondValence=atoms.map(()=>0),bondKeys=new Set();
  atoms.forEach((atom,index)=>{
    if(!atom||typeof atom.element!=='string'||!Object.hasOwn(ATOMIC_WEIGHTS,atom.element))throw new ChemistryError(`L’elemento dell’atomo ${index+1} non è supportato nel disegno.`,'GRAPH_ELEMENT');
    atom.charge=atom.charge??0;
    if(!graphValences(atom.element,atom.charge))throw new ChemistryError(`La carica dell’atomo ${index+1} non è supportata.`,'GRAPH_CHARGE');
    ['x','y','z','x2d','y2d'].forEach(key=>{if(atom[key]!==undefined&&!Number.isFinite(atom[key]))throw new ChemistryError(`La coordinata ${key} dell’atomo ${index+1} non è valida.`,'GRAPH_COORDINATES');});
  });
  bonds.forEach((bond,index)=>{
    if(!bond||!Number.isInteger(bond.from)||!Number.isInteger(bond.to)||bond.from<0||bond.to<0||bond.from>=atoms.length||bond.to>=atoms.length||bond.from===bond.to||![1,2,3].includes(bond.order))throw new ChemistryError(`Il legame ${index+1} deve collegare due atomi diversi con ordine 1, 2 o 3.`,'GRAPH_BOND');
    const key=[bond.from,bond.to].sort((a,b)=>a-b).join(':');
    if(bondKeys.has(key))throw new ChemistryError('Due atomi non possono avere due legami separati: modifica l’ordine del legame esistente.','GRAPH_BOND');
    bondKeys.add(key);bondValence[bond.from]+=bond.order;bondValence[bond.to]+=bond.order;
  });
  const graph=adjacency({atoms,bonds}),visited=new Set([0]),stack=[0];
  while(stack.length){const node=stack.pop();graph[node].forEach(n=>{if(!visited.has(n.to)){visited.add(n.to);stack.push(n.to);}});}
  if(visited.size!==atoms.length)throw new ChemistryError('Collega tutti gli atomi: il disegno deve rappresentare una sola molecola.','GRAPH_DISCONNECTED');
  atoms.forEach((atom,index)=>{
    const allowed=graphValences(atom.element,atom.charge),attached=bondValence[index];
    if(attached>allowed.at(-1))throw new ChemistryError(`L’atomo ${index+1} (${atom.element}) avrebbe troppi legami: ${attached}.`,'GRAPH_VALENCE');
    if(atom.hydrogens===undefined) {
      const target=allowed.find(value=>value>=attached);
      atom.hydrogens=atom.element==='H'?0:target-attached;
      if(!allowed.includes(attached+atom.hydrogens))throw new ChemistryError(`L’atomo ${index+1} (${atom.element}) richiede un legame per questa carica.`,'GRAPH_VALENCE');
    } else {
      if(!Number.isInteger(atom.hydrogens)||atom.hydrogens<0||atom.element==='H'&&atom.hydrogens!==0)throw new ChemistryError(`Gli idrogeni dell’atomo ${index+1} devono essere un intero non negativo.`,'GRAPH_HYDROGENS');
      if(!allowed.includes(attached+atom.hydrogens))throw new ChemistryError(`L’atomo ${index+1} (${atom.element}) ha ${atom.hydrogens} idrogeni e ${attached} legami: questa valenza non è supportata.`,'GRAPH_HYDROGENS');
    }
  });
  let inputBackbone;
  if(graphInput.backbone!==undefined) {
    inputBackbone=graphInput.backbone;
    if(!Array.isArray(inputBackbone)||!inputBackbone.length||new Set(inputBackbone).size!==inputBackbone.length||inputBackbone.some((index,i)=>!Number.isInteger(index)||index<0||index>=atoms.length||i>0&&!graph[inputBackbone[i-1]].some(n=>n.to===index)))throw new ChemistryError('La catena indicata deve essere un percorso continuo di atomi diversi.','GRAPH_BACKBONE');
  }
  const derived=interpretedTreeName({atoms,bonds}),layoutBackbone=inputBackbone||derived?.backbone||longestPath({atoms,bonds});
  const fallbackPositions=atoms.map((atom,index)=>{
    const position=layoutBackbone?.indexOf(index)??-1;
    return position>=0?{x:(position-(layoutBackbone.length-1)/2)*1.3,y:position%2?.39:-.39}:{x:Math.cos(index*2*Math.PI/atoms.length)*Math.max(1.5,atoms.length*.22),y:Math.sin(index*2*Math.PI/atoms.length)*Math.max(1.5,atoms.length*.22)};
  });
  atoms.forEach((atom,index)=>{
    atom.x2d=atom.x2d??atom.x??fallbackPositions[index].x;atom.y2d=atom.y2d??atom.y??fallbackPositions[index].y;
    atom.x=atom.x??atom.x2d;atom.y=atom.y??atom.y2d;atom.z=atom.z??0;
  });
  const totalCharge=atoms.reduce((sum,atom)=>sum+atom.charge,0);
  const formula=formulaFromGraph(atoms)+(totalCharge?(totalCharge>0?'+':'-')+(Math.abs(totalCharge)===1?'':Math.abs(totalCharge)):'');
  const compound={...graphInput,atoms,bonds,name:options.name?cleanName(options.name):derived?.normalizedName||'Molecola disegnata',canonicalName:derived?.canonicalName||'',canonicalNameVerified:false,namingStatus:derived?'reconstructed':'drawn',normalizedName:derived?.normalizedName||'',inputInterpretation:derived?`Struttura ricostruita come ${derived.normalizedName}.`:'Struttura ricostruita dal disegno; nome non determinato.',formula,totalCharge,mass:Number(graphMass(atoms).toFixed(3)),family:inferFamily(atoms,bonds),sourceLabel:options.sourceLabel||'Disegno 2D',conformerSource:'Geometria 3D stimata · disegno non ottimizzato',notes:['La struttura deriva dal disegno; la geometria 3D è stimata e non ottimizzata.'],condensed:null};
  delete compound.sdf3d;delete compound.sourceUrl;delete compound.smiles;delete compound.backbone;
  if(inputBackbone)compound.inputBackbone=[...inputBackbone];
  if(derived)compound.backbone=[...derived.backbone];
  else if(inputBackbone?.every(index=>atoms[index].element==='C'))compound.backbone=[...inputBackbone];
  if(Array.isArray(options.notes))compound.notes.push(...options.notes.filter(note=>typeof note==='string'));
  const condensedBackbone=derived?.backbone||layoutBackbone;
  if(condensedBackbone&&bonds.length===atoms.length-1)compound.condensed=chainFragment(compound,condensedBackbone);
  if(!compound.condensed)compound.notes.push('Per questa struttura ciclica o complessa la condensata non viene ricostruita automaticamente.');
  return compound;
}

function unit(v) {const length=Math.hypot(v.x,v.y,v.z)||1;return{x:v.x/length,y:v.y/length,z:v.z/length};}
function dot(a,b) {return a.x*b.x+a.y*b.y+a.z*b.z;}
function bondDirections(atom,graph,atoms) {return graph.map(n=>unit({x:atoms[n.to].x-atom.x,y:atoms[n.to].y-atom.y,z:atoms[n.to].z-atom.z}));}
const SPHERE_DIRECTIONS=Array.from({length:80},(_,i)=>{const z=1-2*(i+.5)/80,angle=i*Math.PI*(3-Math.sqrt(5));return{x:Math.cos(angle)*Math.sqrt(1-z*z),y:Math.sin(angle)*Math.sqrt(1-z*z),z};});

/** Explicit Lewis graph with H, bond orders and lone pairs; approximate H geometry. */
export function buildLewisGraph(compound) {
  const source=compound.atoms,graph=adjacency(compound);
  const atoms=source.map(a=>({...a,x2d:(a.x2d??a.x)*2.4,y2d:(a.y2d??a.y)*2.4,lonePairs:0}));
  const bonds=compound.bonds.map(b=>({...b}));
  let supported=compound.bonds.every(b=>[1,2,3].includes(b.order));
  source.forEach((atom,index)=>{
    const valence=graph[index].reduce((sum,n)=>sum+n.order,0)+(atom.hydrogens||0);
    const electrons=VALENCE_ELECTRONS[atom.element];
    const lonePairs=(electrons-(atom.charge||0)-valence)/2;
    if(electrons===undefined||lonePairs<0||!Number.isInteger(lonePairs)) {atoms[index].lonePairs=null;supported=false;}
    else atoms[index].lonePairs=lonePairs;
    const takenAngles=graph[index].map(n=>Math.atan2(source[n.to].y2d-atom.y2d,source[n.to].x2d-atom.x2d));
    const directions=bondDirections(atom,graph[index],source);
    for(let h=0;h<(atom.hydrogens||0);h++) {
      const candidates=Array.from({length:24},(_,i)=>i*Math.PI/12);
      const angularDistance=(a,b)=>Math.acos(Math.cos(a-b));
      let angle=candidates[0],score=-Infinity;
      for(const candidate of candidates) {
        let value=takenAngles.length?Math.min(...takenAngles.map(a=>angularDistance(candidate,a))):1;
        // Avoid putting an H label into a heavy atom or an existing H label.
        const px=atoms[index].x2d+Math.cos(candidate)*1.08,py=atoms[index].y2d+Math.sin(candidate)*1.08;
        atoms.forEach((other,otherIndex)=>{if(otherIndex!==index&&Math.hypot(other.x2d-px,other.y2d-py)<.6)value-=3;});
        if(value>score){score=value;angle=candidate;}
      }
      takenAngles.push(angle);
      let vector=SPHERE_DIRECTIONS[0],vectorScore=-Infinity;
      // For saturated atoms target the tetrahedral cosine -1/3; planar sp2 bonds prefer 120 degrees.
      const target=graph[index].some(n=>n.order===2)?-.5:graph[index].some(n=>n.order===3)?-1:-1/3;
      SPHERE_DIRECTIONS.forEach(candidate=>{const value=-directions.reduce((sum,d)=>sum+(dot(candidate,d)-target)**2,0);if(value>vectorScore){vectorScore=value;vector=candidate;}});
      directions.push(vector);
      const hydrogenIndex=atoms.length;
      atoms.push({element:'H',hydrogens:0,charge:0,lonePairs:0,x:atom.x+vector.x*1.09,y:atom.y+vector.y*1.09,z:atom.z+vector.z*1.09,x2d:atoms[index].x2d+Math.cos(angle)*1.08,y2d:atoms[index].y2d+Math.sin(angle)*1.08});
      bonds.push({from:index,to:hydrogenIndex,order:1});
    }
  });
  return {atoms,bonds,supported,notes:supported?[]:['Lewis completa non disponibile per elementi o valenze non supportati.']};
}

const BUNDLED_NAMES = [
  {cid:21201047,names:['5-isobutildecano','5-isobutil-decano','5-isobutyl-decane','5-isobutyldecane','5-(2-metilpropil)decano','5-(2-metilpropil)-decano','5-(2-methylpropyl)decane','5-(2-methylpropyl)-decane'],family:'Alcano ramificato',condensed:'CH3-(CH2)3-CH(CH2-CH(CH3)2)-(CH2)4-CH3'},
  {cid:702,names:['etanolo','ethanol','alcooletilico','alcoletilico'],family:'Alcol',condensed:'CH3-CH2-OH'},
  {cid:241,names:['benzene','benzolo'],family:'Idrocarburo aromatico',condensed:'ciclo-(CH=CH)3'},
  {cid:180,names:['acetone','propanone','propan-2-one','2-propanone'],family:'Chetone',condensed:'CH3-C(=O)-CH3'},
  {cid:176,names:['acidoacetico','acidoetanoico','aceticacid','ethanoicacid'],family:'Acido carbossilico',condensed:'CH3-C(=O)-OH'},
  {cid:2519,names:['caffeina','caffeine','1,3,7-trimethylpurine-2,6-dione'],family:'Alcaloide',condensed:null},
  {cid:962,names:['acqua','water','h2o'],family:'Composto inorganico',condensed:'H-O-H'},
  {cid:887,names:['metanolo','methanol','alcolmetilico','alcoolmetilico'],family:'Alcol',condensed:'CH3-OH'},
];
function resolveBundled(name,key) {
  const definition=BUNDLED_NAMES.find(item=>item.names.includes(key));
  if(!definition)return null;
  const record=BUNDLED_REFERENCE_DATA.find(item=>item.properties.CID===definition.cid);
  const p=record.properties, graph=readPubChem(record.two);
  const compound={name,canonicalName:p.IUPACName,canonicalNameVerified:true,namingStatus:'verified',normalizedName:p.IUPACName,formula:p.MolecularFormula,mass:Number(p.MolecularWeight),family:definition.family,sourceLabel:'PubChem · integrato',sourceUrl:`https://pubchem.ncbi.nlm.nih.gov/compound/${p.CID}`,smiles:p.ConnectivitySMILES,atoms:graph.atoms,bonds:graph.bonds,condensed:definition.condensed,conformerSource:'Conformero 3D PubChem',notes:[]};
  merge3d(compound,readPubChem(record.three));
  const backbone=longestPath(compound);
  if(backbone&&compound.atoms.every(a=>a.element==='C')) {
    const degrees=adjacency(compound);
    const locants=backbone.flatMap((id,i)=>degrees[id].length>2?[i+1]:[]);
    const reverse=locants.map(l=>backbone.length+1-l).sort((a,b)=>a-b);
    if(reverse.join(',')<locants.join(','))backbone.reverse();
    compound.backbone=backbone;
  }
  if(definition.cid===2519)compound.notes.push('Per la caffeina, una struttura a due anelli, la formula semisviluppata conserva i collegamenti che una condensata lineare non descriverebbe.');
  if(definition.cid===241)compound.notes.push('I doppi legami alternati rappresentano una forma di Kekulé: nel benzene gli elettroni π sono delocalizzati.');
  return compound;
}

const BUNDLED_REFERENCE_DATA = [{"properties":{"CID":21201047,"MolecularFormula":"C14H30","MolecularWeight":"198.39","ConnectivitySMILES":"CCCCCC(CCCC)CC(C)C","IUPACName":"5-(2-methylpropyl)decane"},"two":{"id":{"id":{"cid":21201047}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44],"element":[6,6,6,6,6,6,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},"bonds":{"aid1":[1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14],"aid2":[2,3,4,15,5,16,17,6,18,19,7,20,21,9,10,22,8,23,24,11,25,26,12,27,28,29,30,31,32,33,34,13,35,36,14,37,38,39,40,41,42,43,44],"order":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44],"conformers":[{"x":[6.3301,7.1962,5.4641,6.3301,8.0622,4.5981,5.4641,3.732,8.9282,8.0622,5.4641,2.866,4.5981,2,6.8671,7.5947,6.7976,5.8626,5.0656,6.5422,6.9407,8.0622,4.1996,4.9966,5.252,4.8535,4.1306,3.3335,9.2382,9.4651,8.6182,7.4422,8.0622,8.6822,5.6762,6.0747,2.4675,3.2646,4.2881,4.0611,4.9081,2.31,1.4631,1.69],"y":[1.25,1.75,1.75,0.25,1.25,1.25,-0.25,1.75,1.75,0.25,-1.25,1.25,-1.75,1.75,0.94,2.225,2.225,2.225,2.225,-0.3326,0.3577,1.87,0.7751,0.7751,0.3326,-0.3577,2.225,2.225,1.2131,2.06,2.2869,0.25,-0.37,0.25,-1.8326,-1.1423,0.7751,0.7751,-1.2131,-2.06,-2.2869,2.2869,2.06,1.2131]}]}]},"three":{"id":{"id":{"cid":21201047}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44],"element":[6,6,6,6,6,6,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},"bonds":{"aid1":[1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14],"aid2":[2,3,4,15,5,16,17,6,18,19,7,20,21,9,10,22,8,23,24,11,25,26,12,27,28,29,30,31,32,33,34,13,35,36,14,37,38,39,40,41,42,43,44],"order":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44],"conformers":[{"x":[0.5684,1.7598,-0.7802,0.7245,3.1604,-2.0144,0.7615,-3.3032,3.1894,4.2658,0.995,-4.533,1.0233,-5.8173,0.5459,1.7746,1.5683,-0.8581,-0.7857,-0.0693,1.6445,3.3834,-2.0403,-1.9727,-0.1822,1.5585,-3.2638,-3.3854,2.8783,4.2012,2.5392,4.1239,4.2721,5.2504,1.9443,0.201,-4.4668,-4.5689,1.8281,0.0753,1.1905,-5.8276,-6.6838,-5.9303],"y":[-0.7698,-1.1128,-0.9109,0.6524,-1.143,-0.6193,1.8171,-0.7466,-2.1954,-1.4377,3.1771,-0.4242,4.2946,-0.5816,-1.4836,-0.4237,-2.1024,-1.9376,-0.2602,0.8366,0.6817,-0.1641,-1.3087,0.398,1.8532,1.6533,-0.0672,-1.7659,-3.1785,-2.2963,-1.9173,-2.4179,-0.6824,-1.4286,3.1665,3.3771,0.6029,-1.0884,4.1389,4.3505,5.2585,0.0917,-0.3454,-1.608],"z":[-0.2757,0.6555,0.4737,-0.8718,0.0075,-0.3883,0.1271,0.4286,-1.1059,1.0221,-0.5315,-0.4198,0.4999,0.3785,-1.1073,1.5085,1.0942,0.8551,1.3566,-1.6045,-1.4686,-0.4308,-1.2403,-0.7895,0.6829,0.8601,1.2887,0.825,-0.7367,-1.514,-1.9404,1.4897,1.815,0.5427,-1.0792,-1.2602,-0.7971,-1.2912,1.2256,1.0443,0.0092,1.2415,-0.247,0.7417]}]}]}},{"properties":{"CID":702,"MolecularFormula":"C2H6O","MolecularWeight":"46.07","ConnectivitySMILES":"CCO","IUPACName":"ethanol"},"two":{"id":{"id":{"cid":702}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9],"element":[8,6,6,1,1,1,1,1,1]},"bonds":{"aid1":[1,1,2,2,2,3,3,3],"aid2":[2,9,3,4,5,6,7,8],"order":[1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9],"conformers":[{"x":[3.7321,2.866,2,2.4675,3.2646,2.31,1.4631,1.69,4.269],"y":[0.25,-0.25,0.25,-0.7249,-0.7249,0.7869,0.56,-0.2869,-0.06]}]}]},"three":{"id":{"id":{"cid":702}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9],"element":[8,6,6,1,1,1,1,1,1]},"bonds":{"aid1":[1,1,2,2,2,3,3,3],"aid2":[2,9,3,4,5,6,7,8],"order":[1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9],"conformers":[{"x":[-1.1712,-0.0463,1.2175,-0.0958,-0.0952,2.105,1.2426,1.2616,-1.1291],"y":[0.2997,-0.5665,0.2668,-1.212,-1.1938,-0.372,0.9307,0.9052,0.8364],"z":[0,0,0,0.8819,-0.8946,-0.0177,-0.8704,0.8886,0.8099]}]}]}},{"properties":{"CID":241,"MolecularFormula":"C6H6","MolecularWeight":"78.11","ConnectivitySMILES":"C1=CC=CC=C1","IUPACName":"benzene"},"two":{"id":{"id":{"cid":241}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10,11,12],"element":[6,6,6,6,6,6,1,1,1,1,1,1]},"bonds":{"aid1":[1,1,1,2,2,3,3,4,4,5,5,6],"aid2":[2,3,7,4,8,5,9,6,10,6,11,12],"order":[2,1,1,1,1,2,1,2,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10,11,12],"conformers":[{"x":[2.866,3.7321,2,3.7321,2,2.866,2.866,4.269,1.4631,4.269,1.4631,2.866],"y":[1,0.5,0.5,-0.5,-0.5,-1,1.62,0.81,0.81,-0.81,-0.81,-1.62]}]}]},"three":{"id":{"id":{"cid":241}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10,11,12],"element":[6,6,6,6,6,6,1,1,1,1,1,1]},"bonds":{"aid1":[1,1,1,2,2,3,3,4,4,5,5,6],"aid2":[2,3,7,4,8,5,9,6,10,6,11,12],"order":[2,1,1,1,1,2,1,2,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10,11,12],"conformers":[{"x":[-1.2131,-1.2028,-0.0103,0.0104,1.2028,1.2131,-2.1577,-2.1393,-0.0184,0.0184,2.1394,2.1577],"y":[-0.6884,0.7064,-1.3948,1.3948,-0.7063,0.6884,-1.2244,1.2564,-2.4809,2.4808,-1.2563,1.2245],"z":[0,0.0001,0,-0.0001,0,0,0,0.0001,-0.0001,0,0.0001,0]}]}]}},{"properties":{"CID":180,"MolecularFormula":"C3H6O","MolecularWeight":"58.08","ConnectivitySMILES":"CC(=O)C","IUPACName":"propan-2-one"},"two":{"id":{"id":{"cid":180}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10],"element":[8,6,6,6,1,1,1,1,1,1]},"bonds":{"aid1":[1,2,2,3,3,3,4,4,4],"aid2":[2,3,4,5,6,7,8,9,10],"order":[2,1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10],"conformers":[{"x":[3.7321,2.866,2,2.866,2.31,1.4631,1.69,2.246,2.866,3.486],"y":[0.75,0.25,0.75,-0.75,1.2869,1.06,0.2131,-0.75,-1.37,-0.75]}]}]},"three":{"id":{"id":{"cid":180}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10],"element":[8,6,6,6,1,1,1,1,1,1]},"bonds":{"aid1":[1,2,2,3,3,3,4,4,4],"aid2":[2,3,4,5,6,7,8,9,10],"order":[2,1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10],"conformers":[{"x":[0.0003,0,1.281,-1.2813,1.3279,1.326,2.1351,-2.1352,-1.3284,-1.3266],"y":[-1.3171,-0.0872,0.7024,0.7019,1.3235,1.3282,0.0196,0.0187,1.323,1.3278],"z":[-0.0002,0.0006,-0.0002,-0.0002,-0.898,0.8945,0.0027,0.0027,-0.898,0.8945]}]}]}},{"properties":{"CID":176,"MolecularFormula":"C2H4O2","MolecularWeight":"60.05","ConnectivitySMILES":"CC(=O)O","IUPACName":"acetic acid"},"two":{"id":{"id":{"cid":176}},"atoms":{"aid":[1,2,3,4,5,6,7,8],"element":[8,8,6,6,1,1,1,1]},"bonds":{"aid1":[1,1,2,3,3,3,3],"aid2":[4,8,4,4,5,6,7],"order":[1,1,2,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8],"conformers":[{"x":[3.7321,2.866,2,2.866,2.31,1.4631,1.69,4.269],"y":[0.75,-0.75,0.75,0.25,1.2869,1.06,0.2131,0.44]}]}]},"three":{"id":{"id":{"cid":176}},"atoms":{"aid":[1,2,3,4,5,6,7,8],"element":[8,8,6,6,1,1,1,1]},"bonds":{"aid1":[1,1,2,3,3,3,3],"aid2":[4,8,4,4,5,6,7],"order":[1,1,2,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8],"conformers":[{"x":[-0.3035,-0.98,1.3743,-0.0907,1.8368,1.84,1.5207,-1.2598],"y":[1.289,-0.8878,-0.3516,-0.0496,0.057,0.0676,-1.4356,1.5081],"z":[-0.0002,-0.0002,-0.0002,0.0006,-0.9021,0.8952,0.0064,-0.0008]}]}]}},{"properties":{"CID":2519,"MolecularFormula":"C8H10N4O2","MolecularWeight":"194.19","ConnectivitySMILES":"CN1C=NC2=C1C(=O)N(C(=O)N2C)C","IUPACName":"1,3,7-trimethylpurine-2,6-dione"},"two":{"id":{"id":{"cid":2519}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],"element":[8,8,7,7,7,7,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1]},"bonds":{"aid1":[1,2,3,3,3,4,4,4,5,5,5,6,6,7,7,11,12,12,12,13,13,13,14,14,14],"aid2":[9,10,8,10,12,7,11,13,9,10,14,8,11,8,9,15,16,17,18,19,20,21,22,23,24],"order":[2,2,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],"conformers":[{"x":[3.7321,2,3.7321,5.5443,2.866,5.5443,4.5981,4.5981,3.7321,2.866,6.1279,3.7321,5.855,2,6.7479,3.1121,3.7321,4.3521,6.4443,6.0476,5.2656,2.31,1.4631,1.69],"y":[2,-1,-1,0.8047,0.5,-0.8047,0.5,-0.5,1,-0.5,0,-2,1.7553,1,0,-2,-2.62,-2,1.5626,2.3446,1.9479,1.5369,1.31,0.4631]}]}]},"three":{"id":{"id":{"cid":2519}},"atoms":{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],"element":[8,8,7,7,7,7,6,6,6,6,6,6,6,6,1,1,1,1,1,1,1,1,1,1]},"bonds":{"aid1":[1,2,3,3,3,4,4,4,5,5,5,6,6,7,7,11,12,12,12,13,13,13,14,14,14],"aid2":[9,10,8,10,12,7,11,13,9,10,14,8,11,8,9,15,16,17,18,19,20,21,22,23,24],"order":[2,2,1,1,1,1,1,1,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],"conformers":[{"x":[0.47,-3.1271,-0.9686,2.2182,-1.3477,1.4119,0.8579,0.3897,0.0307,-1.9061,2.5032,-1.4276,3.1926,-2.2969,3.5163,-1.0451,-2.5186,-1.0447,4.1992,3.0468,3.0466,-1.8087,-2.9322,-2.9346],"y":[2.5688,-0.4436,-1.3125,0.1412,1.0797,-1.9372,0.2592,-1.0264,1.422,-0.2495,-1.1998,-2.696,1.2061,2.1881,-1.5787,-3.1973,-2.7596,-3.1963,0.7801,1.8092,1.8083,3.1651,2.1027,2.1021],"z":[0.0006,-0.0003,0,-0.0003,-0.0001,0.0002,-0.0008,-0.0004,-0.0006,-0.0004,0.0003,0.0008,0.0003,0.0007,0.0008,-0.8937,0.0011,0.8957,0.0002,-0.8992,0.9004,-0.0003,0.8881,-0.8849]}]}]}},{"properties":{"CID":962,"MolecularFormula":"H2O","MolecularWeight":"18.015","ConnectivitySMILES":"O","IUPACName":"oxidane"},"two":{"id":{"id":{"cid":962}},"atoms":{"aid":[1,2,3],"element":[8,1,1]},"bonds":{"aid1":[1,1],"aid2":[2,3],"order":[1,1]},"coords":[{"aid":[1,2,3],"conformers":[{"x":[2.5369,3.0739,2],"y":[-0.155,0.155,0.155]}]}]},"three":{"id":{"id":{"cid":962}},"atoms":{"aid":[1,2,3],"element":[8,1,1]},"bonds":{"aid1":[1,1],"aid2":[2,3],"order":[1,1]},"coords":[{"aid":[1,2,3],"conformers":[{"x":[0,0.2774,0.6068],"y":[0,0.8929,-0.2383],"z":[0,0.2544,-0.7169]}]}]}},{"properties":{"CID":887,"MolecularFormula":"CH4O","MolecularWeight":"32.042","ConnectivitySMILES":"CO","IUPACName":"methanol"},"two":{"id":{"id":{"cid":887}},"atoms":{"aid":[1,2,3,4,5,6],"element":[8,6,1,1,1,1]},"bonds":{"aid1":[1,1,2,2,2],"aid2":[2,6,3,4,5],"order":[1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6],"conformers":[{"x":[2.5369,3.403,3.713,3.9399,3.093,2],"y":[-0.25,0.25,-0.2869,0.56,0.7869,0.06]}]}]},"three":{"id":{"id":{"cid":887}},"atoms":{"aid":[1,2,3,4,5,6],"element":[8,6,1,1,1,1]},"bonds":{"aid1":[1,1,2,2,2],"aid2":[2,6,3,4,5],"order":[1,1,1,1,1]},"coords":[{"aid":[1,2,3,4,5,6],"conformers":[{"x":[0.7079,-0.7079,-1.0732,-1.0731,-1.0632,0.9936],"y":[0,0,-0.769,-0.1947,0.9786,-0.8804],"z":[0,0,0.6852,-1.0113,0.3312,-0.298]}]}]}}];
