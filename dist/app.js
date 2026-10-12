import { resolveCompound, compoundFromGraph } from './chemistry.js';
import { renderTopology, renderLewis, renderCondensed, create3DViewer, getLewisNotes } from './renderers.js';
import { renderCondensed2D } from './condensed2d.js';
import { createCondensedEditor } from './condensed-editor.js';

const $ = selector => document.querySelector(selector);
const $$ = selector => [...document.querySelectorAll(selector)];
const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const formulaHtml = value => escape(value).replace(/(\d+)/g, '<sub>$1</sub>');
const viewNames = {topology:'Formula topologica',lewis:'Struttura di Lewis',condensed:'Struttura condensata',condensed2d:'Condensata 2D','3d':'Modello 3D'};
const colors = {C:'#4f7466',H:'#d5e1db',O:'#e37482',N:'#789fda',F:'#77bf9c',Cl:'#77bf9c',Br:'#b68b66',I:'#b891c5',S:'#dfbc55',P:'#e5a36b'};
let current, viewer, expandedViewer, request, searchSequence = 0, currentView = 'all', spin = false, toastTimer;
let history = [];
let editor;
try { history = JSON.parse(localStorage.getItem('carbonsketch-history') || '[]').filter(x => typeof x === 'string').slice(0,6); } catch {}

function toast(message) {
  $('#toast').textContent = message;
  $('#toast').hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { $('#toast').hidden = true; }, 3200);
}

function status(message, kind) {
  const node = $('#search-status');
  node.textContent = message;
  node.className = `search-status ${kind || ''}`;
  node.hidden = !message;
}

function renderHistory() {
  const items = history.filter(name => name.toLowerCase() !== current?.name?.toLowerCase());
  $('#history-section').hidden = !items.length;
  $('#history-list').innerHTML = items.map(name => `<button data-compound="${escape(name)}">${escape(name)} ↗</button>`).join('');
}

function presentCompound(compound, {drawing = false} = {}) {
  current = compound;
  $('#compound-name').textContent = compound.name;
  $('#canonical-name').textContent = compound.canonicalName || '';
  $('#canonical-name').hidden = !compound.canonicalName || compound.canonicalName === compound.name;
  $('#molecular-formula').innerHTML = formulaHtml(compound.formula);
  $('#molar-mass').textContent = Number.isFinite(Number(compound.mass)) ? Number(compound.mass).toLocaleString('it-IT',{minimumFractionDigits:2,maximumFractionDigits:2}) : '—';
  $('#molecule-family').textContent = compound.family || 'Composto';
  $('#source-label').textContent = compound.sourceLabel || 'Struttura ricavata dal nome';
  $('#source-link').hidden = !compound.sourceUrl;
  if (compound.sourceUrl && /^https:\/\/(pubchem\.ncbi\.nlm\.nih\.gov|cactus\.nci\.nih\.gov)\//.test(compound.sourceUrl)) $('#source-link').href = compound.sourceUrl;
  else $('#source-link').hidden = true;
  $('#conformer-label').textContent = compound.conformerSource || '3D: geometria illustrativa';
  const notes = [...new Set([...(compound.notes || []),...getLewisNotes(compound)])];
  $('#molecule-notes').textContent = notes.join(' ');
  $('#molecule-notes').hidden = !notes.length;
  $('#name-interpretation').textContent = compound.inputInterpretation || '';
  $('#name-interpretation').hidden = !compound.inputInterpretation;
  updateFormula();
  updateModel();
  status('');
  if (!drawing) {
    history = [compound.name,...history.filter(item=>item.toLowerCase()!==compound.name.toLowerCase())].slice(0,6);
    try { localStorage.setItem('carbonsketch-history', JSON.stringify(history)); } catch {}
  }
  renderHistory();
}

function setInputMode(mode) {
  const drawing = mode === 'draw';
  $('#name-entry').hidden = drawing;
  $('#editor-entry').hidden = !drawing;
  $('#input-mode-name').setAttribute('aria-pressed',String(!drawing));
  $('#input-mode-draw').setAttribute('aria-pressed',String(drawing));
  if (drawing && !editor) {
    editor = createCondensedEditor($('#condensed-editor'), {
      onAnalyze: compound => {
        request?.abort();
        searchSequence++;
        $('#search-button').disabled = false;
        $('#workspace').setAttribute('aria-busy','false');
        presentCompound(compound,{drawing:true});
        setView('condensed2d');
        const url = new URL(location.href);
        url.searchParams.delete('q');
        url.searchParams.set('drawing','1');
        url.searchParams.set('view','condensed2d');
        window.history.replaceState(null,'',url);
        try { localStorage.setItem('carbonsketch-drawing',JSON.stringify({atoms:compound.atoms,bonds:compound.bonds,backbone:compound.backbone})); } catch {}
        toast('Disegno analizzato. Esplora tutte le viste della molecola.');
      },
      onError: message => toast(message)
    });
  }
}

function updateFormula() {
  const diagramSize = selector => {const node=$(selector);return {width:node.clientWidth ? Math.max(280,node.clientWidth-30) : 760,height:node.clientHeight ? Math.max(180,node.clientHeight-20) : 380};};
  $('#numbering-toggle').disabled = !current.backbone?.length;
  $('#numbering-toggle').closest('label').title = current.backbone?.length ? 'Numera i carboni della catena principale' : 'Numerazione della catena principale non disponibile per questa struttura';
  if (!current.backbone?.length) $('#numbering-toggle').checked = false;
  $('#topology-canvas').innerHTML = renderTopology(current, {...diagramSize('#topology-canvas'),numbering:$('#numbering-toggle').checked});
  $('#lewis-canvas').innerHTML = renderLewis(current,diagramSize('#lewis-canvas'));
  const partialLewis = getLewisNotes(current).some(note=>/non disponibile/i.test(note));
  $('.lewis-card h3').textContent = partialLewis ? 'Lewis · vista parziale' : 'Struttura di Lewis';
  $('.condensed-content').classList.toggle('graph-condensed',!current.condensed);
  $('#condensed-formula').innerHTML = current.condensed ? formulaHtml(current.condensed) : renderCondensed(current,diagramSize('#condensed-formula'));
  $('.formula-kicker').textContent = current.condensed ? 'I LEGAMI, IN UNA RIGA' : 'FORMULA SEMISVILUPPATA';
  $('#condensed-description').textContent = current.condensedNote || (current.condensed ? 'Le parentesi indicano gruppi collegati all’atomo che le precede.' : 'Gruppi di atomi espliciti, con tutti i collegamenti della struttura.');
  $('#copy-condensed').disabled = !current.condensed;
  $('#copy-condensed').title = current.condensed ? 'Copia struttura condensata' : 'Per questa struttura è disponibile il download SVG';
  const canvas2d = $('#condensed-2d-canvas');
  const diagram2d = renderCondensed2D(current,{width:Math.max(320,canvas2d.clientWidth-34),height:240,numbering:$('#condensed-2d-numbering').checked});
  canvas2d.innerHTML = diagram2d;
  const svg2d = canvas2d.querySelector('svg');
  svg2d.style.width = `${svg2d.getAttribute('width')}px`;
  svg2d.style.height = `${svg2d.getAttribute('height')}px`;
  $('#condensed-2d-note').textContent = svg2d.dataset.fallback === 'true' ? 'Gruppi di atomi espliciti, con tutti i collegamenti della struttura.' : 'La catena principale è orizzontale; i rami sono collegati al loro carbonio.';
  $('#condensed-2d-numbering').disabled = svg2d.dataset.fallback === 'true';
}

function updateModel() {
  viewer?.destroy();
  viewer = undefined;
  if (currentView !== 'all' && currentView !== '3d') return;
  viewer = create3DViewer($('#model-canvas'), current, {style:$('#model-style').value,hydrogens:$('#hydrogens-toggle').checked});
  viewer.setSpin(spin);
  if (!viewer.available) {
    $('#model-canvas').innerHTML = '<div class="model-fallback"><svg class="icon"><use href="#i-cube"/></svg><strong>Il modello 3D non è disponibile</strong><p>Il browser non ha attivato WebGL. Prova ad abilitare l’accelerazione grafica o a usare un altro browser.</p></div>';
  }
  $('.model-hint').hidden = !viewer.available;
  $('#model-legend').innerHTML = [...new Set(current.atoms.map(a=>a.element).concat(current.atoms.some(a=>a.hydrogens>0)?['H']:[]))].map(element => `<span><i style="background:${colors[element] || '#888'}"></i>${escape(element)}</span>`).join('');
}

async function search(name, {initial = false} = {}) {
  const trimmed = name.trim();
  if (!trimmed) { $('#compound-input').focus(); return; }
  if (trimmed.length > 160) { status('Usa un nome composto da meno di 160 caratteri.', 'error'); return; }
  request?.abort();
  request = new AbortController();
  const sequence = ++searchSequence;
  $('#compound-input').value = trimmed;
  $('#search-button').disabled = true;
  $('#workspace').setAttribute('aria-busy','true');
  status('Sto cercando la struttura del composto…', 'loading');
  try {
    const compound = await resolveCompound(trimmed, {signal:request.signal});
    if (sequence !== searchSequence) return;
    presentCompound(compound);
    const url = new URL(location.href);
    url.searchParams.set('q',trimmed);
    url.searchParams.delete('drawing');
    window.history.replaceState(null,'',url);
    if (!initial && matchMedia('(max-width:760px)').matches) $('#compound-name').scrollIntoView({block:'start',behavior:matchMedia('(prefers-reduced-motion:reduce)').matches?'auto':'smooth'});
  } catch (error) {
    if (error.name === 'AbortError' || sequence !== searchSequence) return;
    status(error.message || 'Non ho trovato il composto. Controlla il nome oppure prova il nome inglese.', 'error');
    if (!current && !initial) $('#compound-name').textContent = 'Nessun composto selezionato';
  } finally {
    if (sequence === searchSequence) {
      $('#search-button').disabled = false;
      $('#workspace').setAttribute('aria-busy','false');
    }
  }
}

function setView(view) {
  currentView = view;
  $$('.view-tab').forEach(tab => {const active = tab.dataset.view === view;tab.classList.toggle('active',active);tab.setAttribute('aria-selected',String(active));tab.tabIndex = active ? 0 : -1;});
  $$('.formula-card').forEach(card => {card.hidden = view !== 'all' && card.dataset.card !== view;});
  $('#formula-grid').classList.toggle('focused',view !== 'all');
  $('#view-count').textContent = view === 'all' ? '5 prospettive, 1 molecola' : viewNames[view];
  const url = new URL(location.href);
  url.searchParams.set('view',view);
  window.history.replaceState(null,'',url);
  if (current) updateFormula();
  if (current && (view === '3d' || view === 'all')) requestAnimationFrame(() => { if (!viewer) updateModel(); else viewer.resize(); });
  if (view !== '3d' && view !== 'all') viewer?.setSpin(false);
  else viewer?.setSpin(spin);
}

function openGuide() { $('#guide-dialog').showModal(); }
function caption(mode) {
  if (mode === '3d') return `${current.conformerSource || 'Geometria illustrativa'} · Trascina per ruotare, usa la rotella o due dita per lo zoom.`;
  if (mode === 'topology') return 'Ogni vertice e ogni estremità senza etichetta rappresenta un atomo di carbonio. Gli idrogeni legati al carbonio sono sottintesi.';
  if (mode === 'lewis') return ['Ogni linea rappresenta una coppia elettronica condivisa; i puntini indicano coppie libere.',...getLewisNotes(current)].join(' ');
  if (mode === 'condensed2d') return $('#condensed-2d-canvas svg')?.dataset.fallback !== 'true' ? 'Catena principale in orizzontale e ramificazioni collegate ai rispettivi carboni. Gli idrogeni sono raggruppati in CH₃, CH₂ e CH. Scorri orizzontalmente per leggere le catene lunghe.' : 'Formula semisviluppata: gli idrogeni sono raggruppati con l’atomo a cui sono legati e tutti i collegamenti della struttura restano espliciti.';
  return current.condensedNote || (current.condensed ? 'Le parentesi mostrano gruppi collegati all’atomo che le precede.' : 'Formula semisviluppata: gli idrogeni sono raggruppati con l’atomo a cui sono legati. Tutti i collegamenti, compresi gli anelli, sono conservati.');
}

function expand(mode) {
  if (!current) return;
  expandedViewer?.destroy();
  expandedViewer = undefined;
  $('#expanded-title').textContent = viewNames[mode];
  $('#expanded-compound').textContent = current.name;
  $('#expanded-caption').textContent = caption(mode);
  $('#expanded-content').replaceChildren();
  $('#expanded-content').classList.toggle('expanded-condensed-2d',mode === 'condensed2d');
  $('#expand-dialog').showModal();
  if (mode === '3d') {
    const container = document.createElement('div');
    container.className = 'expanded-model';
    $('#expanded-content').append(container);
    expandedViewer = create3DViewer(container,current,{style:$('#model-style').value,hydrogens:$('#hydrogens-toggle').checked});
    if (!expandedViewer.available) container.innerHTML = '<div class="model-fallback">Le funzionalità 3D richiedono WebGL.</div>';
  } else if (mode === 'condensed2d') {
    $('#expanded-content').innerHTML = renderCondensed2D(current,{width:Math.max(800,$('#expanded-content').clientWidth-30),height:Math.min(540,$('#expanded-content').clientHeight-20),numbering:$('#condensed-2d-numbering').checked});
    const svg = $('#expanded-content svg');
    svg.style.width = `${svg.getAttribute('width')}px`;
    svg.style.height = `${svg.getAttribute('height')}px`;
  } else if (mode === 'condensed') {
    $('#expanded-content').innerHTML = current.condensed ? `<p class="condensed-formula">${formulaHtml(current.condensed)}</p>` : renderCondensed(current,{width:1100,height:650});
  } else $('#expanded-content').innerHTML = mode === 'lewis' ? renderLewis(current,{width:1100,height:650}) : renderTopology(current,{numbering:$('#numbering-toggle').checked,width:1100,height:650});
}

function saveFile(content, mime, extension, label) {
  const blob = content instanceof Blob ? content : new Blob([content],{type:mime});
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = `${current.name.replace(/[^\p{L}\p{N}-]/gu,'_')}-${label || (extension.includes('svg')?'struttura':'carbonsketch')}.${extension}`;
  document.body.append(link);link.click();link.remove();
  setTimeout(()=>URL.revokeObjectURL(url),1000);
}

function download(mode) {
  if (!current) return;
  if (mode === 'condensed') {
    if (!current.condensed) saveFile(renderCondensed(current,{width:1500,height:850}),'image/svg+xml','svg');
    else saveFile(`${current.name}\n${current.canonicalName || ''}\nFormula molecolare: ${current.formula}\nMassa molare: ${current.mass} g/mol\nStruttura condensata: ${current.condensed}\n${caption(mode)}\nFonte: ${current.sourceUrl || current.sourceLabel}\n`,'text/plain; charset=utf-8','txt');
  } else if (mode === 'condensed2d') {
    saveFile(renderCondensed2D(current,{width:1100,height:380,numbering:$('#condensed-2d-numbering').checked}),'image/svg+xml','svg','condensata-2D');
  } else if (mode === '3d') {
    const uri = viewer?.exportPng();
    if (!uri) { toast('Il modello 3D non è disponibile.'); return; }
    const link = document.createElement('a');link.href = uri;link.download = `${current.name}-3D.png`;document.body.append(link);link.click();link.remove();
  } else {
    const svg = mode === 'lewis' ? renderLewis(current,{width:1500,height:850}) : renderTopology(current,{numbering:$('#numbering-toggle').checked,width:1500,height:850});
    saveFile(svg,'image/svg+xml','svg');
  }
  toast('La struttura è pronta per il download.');
}

$('#search-form').addEventListener('submit',event=>{event.preventDefault();search($('#compound-input').value);});
$('#input-mode-name').addEventListener('click',()=>setInputMode('name'));
$('#input-mode-draw').addEventListener('click',()=>setInputMode('draw'));
$('#editor-import').addEventListener('click',()=>{if(current){setInputMode('draw');editor.loadCompound(current);}});
document.addEventListener('click',event=>{
  const example = event.target.closest('[data-compound]');if(example) {setInputMode('name');search(example.dataset.compound);}
  const tab = event.target.closest('[data-view]');if(tab) setView(tab.dataset.view);
  const expandButton = event.target.closest('[data-expand]');if(expandButton) expand(expandButton.dataset.expand);
  const downloadButton = event.target.closest('[data-download]');if(downloadButton) download(downloadButton.dataset.download);
  const close = event.target.closest('[data-close-dialog]');if(close) close.closest('dialog').close();
});
$('#guide-button').addEventListener('click',openGuide);
$('#learn-button').addEventListener('click',openGuide);
$('#numbering-toggle').addEventListener('change',()=>{if(current) updateFormula();});
$('#condensed-2d-numbering').addEventListener('change',()=>{if(current) updateFormula();});
$('#model-style').addEventListener('change',()=>viewer?.setStyle($('#model-style').value));
$('#hydrogens-toggle').addEventListener('change',()=>viewer?.setHydrogens($('#hydrogens-toggle').checked));
$('#reset-3d').addEventListener('click',()=>{viewer?.reset();toast('Orientamento iniziale ripristinato.');});
$('#spin-button').addEventListener('click',()=>{
  spin = !spin;viewer?.setSpin(spin);
  $('#spin-button').setAttribute('aria-pressed',String(spin));
  $('#spin-button').setAttribute('aria-label',spin?'Ferma rotazione automatica':'Attiva rotazione automatica');
  $('#spin-button use').setAttribute('href',spin?'#i-pause':'#i-play');
});
$('#copy-condensed').addEventListener('click',async()=>{
  if(!current?.condensed){toast('La struttura condensata non è disponibile.');return;}
  try {await navigator.clipboard.writeText(current.condensed);toast('Struttura condensata copiata.');}
  catch {const text=document.createElement('textarea');text.value=current.condensed;text.style.cssText='position:fixed;top:-1000px';document.body.append(text);text.select();const copied=document.execCommand('copy');text.remove();toast(copied?'Struttura condensata copiata.':'Seleziona la formula per copiarla.');}
});
$('#expand-dialog').addEventListener('close',()=>{expandedViewer?.destroy();expandedViewer=undefined;});
$$('dialog').forEach(dialog=>dialog.addEventListener('click',event=>{if(event.target===dialog){const r=dialog.getBoundingClientRect();if(event.clientX<r.left||event.clientX>r.right||event.clientY<r.top||event.clientY>r.bottom)dialog.close();}}));
$('.view-tabs').addEventListener('keydown',event=>{
  if(!['ArrowLeft','ArrowRight','Home','End'].includes(event.key))return;
  event.preventDefault();const tabs=$$('.view-tab');let index=tabs.findIndex(tab=>tab.dataset.view===currentView);index=event.key==='Home'?0:event.key==='End'?tabs.length-1:(index+(event.key==='ArrowRight'?1:-1)+tabs.length)%tabs.length;setView(tabs[index].dataset.view);tabs[index].focus();
});
document.addEventListener('keydown',event=>{if(event.key==='/'&&!event.ctrlKey&&!event.metaKey&&!['INPUT','TEXTAREA','SELECT'].includes(document.activeElement.tagName)&&!document.querySelector('dialog[open]')){event.preventDefault();setInputMode('name');$('#compound-input').focus();$('#compound-input').select();}});
document.addEventListener('visibilitychange',()=>viewer?.setSpin(!document.hidden&&spin));
let resizeFrame;
window.addEventListener('resize',()=>{cancelAnimationFrame(resizeFrame);resizeFrame=requestAnimationFrame(()=>{if(current)updateFormula();});});
renderHistory();
const initialView = new URL(location.href).searchParams.get('view');
setView(Object.hasOwn(viewNames,initialView) ? initialView : 'all');
const initialUrl = new URL(location.href);
let restoredDrawing = false;
if (initialUrl.searchParams.get('drawing') === '1') {
  try {
    const saved = JSON.parse(localStorage.getItem('carbonsketch-drawing'));
    if (saved) {
      presentCompound(compoundFromGraph(saved),{drawing:true});
      setInputMode('draw');
      editor.loadCompound(current);
      restoredDrawing = true;
    }
  } catch {}
}
if (!restoredDrawing) {
  const initialQuery = initialUrl.searchParams.get('q') || '5-isobutil decano';
  search(initialQuery,{initial:true}).then(async()=>{if(searchSequence===1&&!current&&initialQuery!=='5-isobutil decano'){const error=$('#search-status').textContent;await search('5-isobutil decano',{initial:true});status(error+' È mostrato il composto di esempio.','error');}});
}
