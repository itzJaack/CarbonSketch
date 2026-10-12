import { compoundFromGraph } from './chemistry.js';

const MAX_ATOMS = 100;
const GROUPS = ['CH3', 'CH2', 'CH', 'C', 'CH4', 'OH', 'O', 'NH2', 'NH', 'N', 'NH3', 'F', 'Cl', 'Br', 'I', 'SH', 'S'];
const ELEMENTS = ['Cl', 'Br', 'C', 'N', 'O', 'F', 'I', 'S', 'P'];
const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const groupLabel = atom => `${atom.element}${atom.hydrogens ? `H${atom.hydrogens > 1 ? atom.hydrogens : ''}` : ''}`;
const labelMarkup = label => escape(label).replace(/(\d+)/g, '<tspan baseline-shift="sub" font-size="13">$1</tspan>');
const round = value => Number(value.toFixed(2));

function formulaError(message) {
  const error = new Error(message); error.name = 'CondensedFormulaError'; return error;
}

function appendFragment(target, source) {
  if (target.atoms.length + source.atoms.length > MAX_ATOMS) throw formulaError(`Usa al massimo ${MAX_ATOMS} gruppi nella formula.`);
  const offset = target.atoms.length;
  target.atoms.push(...source.atoms.map(atom => ({ ...atom })));
  target.bonds.push(...source.bonds.map(bond => ({ ...bond, from: bond.from + offset, to: bond.to + offset })));
  return source.chain.map(index => index + offset);
}

function placeFormulaGraph(graph) {
  const neighbors = graph.atoms.map(() => []), seen = new Set(graph.backbone);
  graph.bonds.forEach(bond => { neighbors[bond.from].push(bond.to); neighbors[bond.to].push(bond.from); });
  graph.backbone.forEach((index, i) => Object.assign(graph.atoms[index], { x2d: i * 1.4, y2d: 0, x: i * 1.4, y: 0, z: 0 }));
  function visit(index, parent, direction) {
    if (seen.has(index)) return;
    seen.add(index);
    const anchor = graph.atoms[parent];
    Object.assign(graph.atoms[index], { x2d: anchor.x2d + Math.cos(direction) * 1.4, y2d: anchor.y2d + Math.sin(direction) * 1.4, z: 0 });
    graph.atoms[index].x = graph.atoms[index].x2d; graph.atoms[index].y = graph.atoms[index].y2d;
    const children = neighbors[index].filter(next => next !== parent && !seen.has(next));
    children.forEach((next, i) => visit(next, index, direction + (i - (children.length - 1) / 2) * .8));
  }
  graph.backbone.forEach(index => neighbors[index].filter(next => !seen.has(next)).forEach((next, i) => visit(next, index, i % 2 ? -Math.PI / 2 : Math.PI / 2)));
  return graph;
}

/** Parse an explicit condensed formula. Written H counts are checked, never silently changed. */
export function parseCondensedFormula(value) {
  if (typeof value !== 'string' || !value.trim()) throw formulaError('Scrivi una formula condensata, per esempio CH3-CH2-CH3.');
  if (value.length > 1800) throw formulaError('La formula è troppo lunga. Usa al massimo 100 gruppi.');
  const text = value.trim().replace(/[₀₁₂₃₄₅₆₇₈₉]/g, character => '₀₁₂₃₄₅₆₇₈₉'.indexOf(character)).replace(/[‐‑–—−]/g, '-').replace(/[≡#]/g, '#').replace(/\s+/g, '');
  let cursor = 0;
  const readCount = () => {
    const match = text.slice(cursor).match(/^\d+/);
    if (!match) return 1;
    cursor += match[0].length;
    const count = Number(match[0]);
    if (!Number.isInteger(count) || count < 1 || count > MAX_ATOMS) throw formulaError('Il numero di ripetizioni deve essere compreso tra 1 e 100.');
    return count;
  };
  function readSequence(depth = 0) {
    if (depth > 20) throw formulaError('Usa al massimo 20 livelli di parentesi.');
    const fragment = { atoms: [], bonds: [], chain: [] };
    let previous = null, pendingOrder = 1, explicitBond = false;
    while (cursor < text.length && text[cursor] !== ')') {
      const character = text[cursor];
      if (['-', '=', '#'].includes(character)) {
        if (previous === null || explicitBond) throw formulaError('Ogni simbolo di legame deve stare fra due gruppi.');
        pendingOrder = character === '=' ? 2 : character === '#' ? 3 : 1;
        explicitBond = true; cursor++; continue;
      }
      if (character === '(') {
        cursor++;
        let attachmentOrder = pendingOrder;
        if (['-', '=', '#'].includes(text[cursor])) {
          if (previous === null || explicitBond) throw formulaError('Il legame di un gruppo fra parentesi deve seguire un atomo e non un altro legame.');
          attachmentOrder = text[cursor] === '=' ? 2 : text[cursor] === '#' ? 3 : 1;
          cursor++;
        }
        const inner = readSequence(depth + 1);
        if (text[cursor] !== ')') throw formulaError('Manca una parentesi chiusa nella formula.');
        cursor++;
        if (!inner.atoms.length) throw formulaError('Le parentesi non possono essere vuote.');
        const count = readCount();
        // A repeated methylene segment is a chain, while C(CH3)3 is three branches.
        const repeatedMethylene = inner.atoms.length === 1 && inner.atoms[0].element === 'C' && inner.atoms[0].hydrogens === 2;
        const isChain = previous === null || explicitBond || repeatedMethylene;
        for (let copy = 0; copy < count; copy++) {
          const mapped = appendFragment(fragment, inner), first = mapped[0], last = mapped.at(-1);
          if (previous !== null) fragment.bonds.push({ from: previous, to: first, order: attachmentOrder });
          if (isChain) { fragment.chain.push(...mapped); previous = last; }
        }
        pendingOrder = 1; explicitBond = false; continue;
      }
      const element = ELEMENTS.find(symbol => text.startsWith(symbol, cursor));
      if (!element) throw formulaError(`Gruppo non riconosciuto vicino a «${text.slice(cursor, cursor + 12)}». Usa gruppi come CH3, CH2, CH, C, OH o O.`);
      cursor += element.length;
      let hydrogens = 0;
      if (text[cursor] === 'H') {
        cursor++;
        const match = text.slice(cursor).match(/^\d+/);
        if (match) { cursor += match[0].length; hydrogens = Number(match[0]); }
        else hydrogens = 1;
        if (hydrogens < 1 || hydrogens > 4 || !Number.isInteger(hydrogens)) throw formulaError('Ogni gruppo può dichiarare da 1 a 4 idrogeni; senza H il numero è zero.');
      }
      if (fragment.atoms.length >= MAX_ATOMS) throw formulaError(`Usa al massimo ${MAX_ATOMS} gruppi nella formula.`);
      const index = fragment.atoms.length;
      fragment.atoms.push({ element, hydrogens, charge: 0, group: `${element}${hydrogens ? `H${hydrogens > 1 ? hydrogens : ''}` : ''}` });
      fragment.chain.push(index);
      if (previous !== null) fragment.bonds.push({ from: previous, to: index, order: pendingOrder });
      previous = index; pendingOrder = 1; explicitBond = false;
    }
    if (explicitBond) throw formulaError('La formula termina con un legame senza gruppo.');
    return fragment;
  }
  const parsed = readSequence();
  if (cursor !== text.length) throw formulaError('C’è una parentesi chiusa senza apertura nella formula.');
  const graph = placeFormulaGraph({ atoms: parsed.atoms, bonds: parsed.bonds, backbone: parsed.chain });
  compoundFromGraph(graph); // Shared scientific validation checks every explicitly written H count.
  return graph;
}

function bondSvg(bond, nodes) {
  const from = nodes.find(node => node.id === bond.from), to = nodes.find(node => node.id === bond.to);
  if (!from || !to) return '';
  const dx = to.x - from.x, dy = to.y - from.y, length = Math.hypot(dx, dy) || 1;
  const ux = dx / length, uy = dy / length, inset = Math.min(35, length * .26);
  const offsets = bond.order === 3 ? [-5, 0, 5] : bond.order === 2 ? [-3, 3] : [0];
  const visible = offsets.map(offset => `<line x1="${round(from.x + ux * inset - uy * offset)}" y1="${round(from.y + uy * inset + ux * offset)}" x2="${round(to.x - ux * inset - uy * offset)}" y2="${round(to.y - uy * inset + ux * offset)}"/>`).join('');
  return `<g class="ce-bond" data-ce-bond="${bond.from}-${bond.to}" role="button" tabindex="0" aria-label="Legame ${bond.order}, gruppo ${bond.from} con ${bond.to}. In modalità Elimina premi per rimuoverlo.">${visible}<line class="ce-bond-hit" x1="${round(from.x)}" y1="${round(from.y)}" x2="${round(to.x)}" y2="${round(to.y)}" stroke="transparent" stroke-width="18"/></g>`;
}

/** A self-contained accessible group-and-bond drawing editor. */
export function createCondensedEditor(container, { onAnalyze = () => {}, onError = () => {} } = {}) {
  if (!container) throw new Error('Il contenitore del disegno non è disponibile.');
  let nodes = [], bonds = [], selected = null, linkStart = null, mode = 'add', nextId = 1, destroyed = false, backbone = [];
  let width = 960, height = 420;
  container.dataset.condensedEditor = '';
  container.innerHTML = `<div class="ce-toolbar"><label class="ce-control">Gruppo<select data-ce-group aria-label="Gruppo da disegnare">${GROUPS.map(group => `<option value="${group}">${group}</option>`).join('')}</select></label><label class="ce-control">Legame<select data-ce-order aria-label="Ordine del legame"><option value="1">Singolo —</option><option value="2">Doppio =</option><option value="3">Triplo ≡</option></select></label><div class="ce-modes" role="group" aria-label="Strumenti del disegno"><button type="button" class="ce-mode active" data-ce-mode="add" aria-pressed="true">Aggiungi</button><button type="button" class="ce-mode" data-ce-mode="link" aria-pressed="false">Collega</button><button type="button" class="ce-mode" data-ce-mode="move" aria-pressed="false">Sposta</button><button type="button" class="ce-mode" data-ce-mode="delete" aria-pressed="false">Elimina</button></div><button type="button" data-ce-add>Aggiungi gruppo</button><button type="button" data-ce-reset>Nuovo disegno</button></div><p class="ce-instructions">Scegli un gruppo e tocca il foglio. Con <strong>Collega</strong> tocca due gruppi; con <strong>Sposta</strong> scegli un gruppo e poi la nuova posizione. Puoi anche usare i controlli sotto il foglio.</p><div class="ce-stage"><svg class="ce-svg" data-ce-canvas xmlns="http://www.w3.org/2000/svg" width="960" height="420" viewBox="0 0 960 420" role="group" aria-label="Foglio per disegnare la formula condensata 2D"></svg></div><div class="ce-status" data-ce-status role="status" aria-live="polite">Il foglio è vuoto. Aggiungi il primo gruppo.</div><div class="ce-properties"><label class="ce-control">Primo gruppo<select data-ce-from aria-label="Seleziona primo gruppo"></select></label><label class="ce-control">Secondo gruppo<select data-ce-to aria-label="Seleziona secondo gruppo"></select></label><button type="button" data-ce-connect>Collega i gruppi</button><button type="button" data-ce-change>Modifica il primo gruppo</button><button type="button" data-ce-delete>Elimina il primo gruppo</button></div><div class="ce-formula-row"><label class="ce-control">Oppure scrivi la condensata<textarea data-ce-formula rows="2" spellcheck="false" aria-label="Formula condensata da disegnare" placeholder="CH3-CH2-CH(CH2-CH3)-CH2-CH(C(CH3)3)-CH2-CH2-CH3"></textarea></label><button type="button" data-ce-read>Leggi formula</button></div><div class="ce-actions"><p>Gli idrogeni restano quelli scritti: se cambi i legami, modifica anche il gruppo (CH3, CH2, CH o C).</p><button type="button" class="ce-analyze" data-ce-analyze>Analizza disegno</button></div>`;
  const $ = selector => container.querySelector(selector);
  const canvas = $('[data-ce-canvas]'), status = $('[data-ce-status]');
  const notify = (message, error = false) => { status.textContent = message; status.classList.toggle('error', error); if (error) onError(message); };
  const groupAtom = () => {
    const label = $('[data-ce-group]').value, match = label.match(/^(Cl|Br|[CNOFISP])(?:H(\d*))?$/);
    return { element: match[1], hydrogens: match[2] !== undefined ? Number(match[2] || 1) : 0, charge: 0, label };
  };
  function render() {
    if (destroyed) return;
    width = Math.max(960, ...nodes.map(node => node.x + 54)); height = Math.max(420, ...nodes.map(node => node.y + 54));
    canvas.setAttribute('viewBox', `0 0 ${width} ${height}`); canvas.setAttribute('width', width); canvas.setAttribute('height', height);
    canvas.style.width = `${width}px`; canvas.style.height = `${height}px`;
    canvas.innerHTML = `<defs><pattern id="ce-grid" width="24" height="24" patternUnits="userSpaceOnUse"><circle cx="1" cy="1" r="1" fill="#d8e3dc"/></pattern></defs><rect width="100%" height="100%" fill="url(#ce-grid)"/><g class="ce-bonds">${bonds.map(bond => bondSvg(bond, nodes)).join('')}</g><g class="ce-nodes">${nodes.map(node => `<g class="ce-node${selected === node.id || linkStart === node.id ? ' selected' : ''}" data-ce-node="${node.id}" transform="translate(${round(node.x)} ${round(node.y)})" role="button" tabindex="0" aria-pressed="${selected === node.id}" aria-label="Gruppo ${node.id}: ${escape(node.label)}"><rect x="-35" y="-24" width="70" height="48" rx="12"/><text x="0" y="1" text-anchor="middle" dominant-baseline="central" font-size="21">${labelMarkup(node.label)}</text><text class="ce-node-index" x="0" y="-32" text-anchor="middle" font-size="11">${node.id}</text></g>`).join('')}</g>${!nodes.length ? '<text class="ce-empty" x="480" y="210" text-anchor="middle">Tocca qui per aggiungere un gruppo</text>' : ''}`;
    for (const selector of ['[data-ce-from]', '[data-ce-to]']) {
      const select = $(selector), previous = select.value;
      select.innerHTML = `<option value="">Scegli un gruppo</option>${nodes.map(node => `<option value="${node.id}">${node.id} · ${escape(node.label)}</option>`).join('')}`;
      select.value = nodes.some(node => String(node.id) === previous) ? previous : '';
    }
    if (selected !== null) $('[data-ce-from]').value = String(selected);
    $('[data-ce-delete]').disabled = !nodes.length; $('[data-ce-change]').disabled = !nodes.length;
    $('[data-ce-connect]').disabled = nodes.length < 2; $('[data-ce-analyze]').disabled = !nodes.length;
  }
  function add(x, y) {
    if (nodes.length >= MAX_ATOMS) { notify('Usa al massimo 100 gruppi nel disegno.', true); return; }
    const node = { id: nextId++, ...groupAtom(), x: Math.max(42, x), y: Math.max(46, y) };
    nodes.push(node); selected = node.id; backbone = []; render(); notify(`Aggiunto ${node.label}. Usa Collega per unirlo agli altri gruppi.`);
  }
  function remove(id) {
    nodes = nodes.filter(node => node.id !== id); bonds = bonds.filter(bond => bond.from !== id && bond.to !== id);
    selected = null; linkStart = null; backbone = []; render(); notify('Gruppo eliminato insieme ai suoi legami.');
  }
  function connect(from, to) {
    if (!from || !to || from === to || !nodes.some(node => node.id === from) || !nodes.some(node => node.id === to)) { notify('Scegli due gruppi diversi da collegare.', true); return; }
    const order = Number($('[data-ce-order]').value), existing = bonds.find(bond => bond.from === from && bond.to === to || bond.from === to && bond.to === from);
    if (existing) existing.order = order;
    else bonds.push({ from, to, order });
    backbone = []; linkStart = null; selected = to; render(); notify(existing ? 'Ordine del legame aggiornato. Controlla anche gli idrogeni dei due gruppi.' : 'Legame aggiunto. Controlla gli idrogeni e premi Analizza disegno.');
  }
  function useNode(id) {
    if (mode === 'delete') { remove(id); return; }
    if (mode === 'link') {
      if (linkStart !== null) { connect(linkStart, id); return; }
      linkStart = id; selected = id; render(); notify('Primo gruppo selezionato. Tocca il secondo gruppo.'); return;
    }
    selected = id; render(); notify(mode === 'move' ? 'Tocca il foglio nel punto in cui vuoi spostare questo gruppo.' : `Selezionato il gruppo ${id}.`);
  }
  function drawingGraph() {
    const indices = new Map(nodes.map((node, index) => [node.id, index]));
    return { atoms: nodes.map(node => ({ element: node.element, hydrogens: node.hydrogens, charge: node.charge || 0, group: node.label, x2d: (node.x - 480) / 70, y2d: -(node.y - 210) / 70, x: (node.x - 480) / 70, y: -(node.y - 210) / 70, z: 0 })), bonds: bonds.map(bond => ({ from: indices.get(bond.from), to: indices.get(bond.to), order: bond.order })), ...(backbone.length ? { backbone: backbone.map(id => indices.get(id)) } : {}) };
  }
  function loadCompound(compound) {
    if (destroyed) return;
    const atoms = compound.atoms || [];
    if (!atoms.length || atoms.length > MAX_ATOMS) { notify('Questa struttura non può essere caricata: usa da 1 a 100 gruppi.', true); return; }
    if (atoms.some(atom => !ELEMENTS.includes(atom.element) || atom.element === 'H' || atom.charge)) { notify('Il disegno supporta C, N, O, F, Cl, Br, I, S e P neutri. Questa struttura contiene gruppi non supportati.', true); return; }
    const points = atoms.map((atom, index) => ({ x: Number.isFinite(Number(atom.x2d)) ? Number(atom.x2d) : index * 1.4, y: Number.isFinite(Number(atom.y2d)) ? -Number(atom.y2d) : 0 }));
    // Keep even extreme imported coordinates finite before calculating spans.
    const magnitude = Math.max(1, ...points.flatMap(point => [Math.abs(point.x), Math.abs(point.y)]));
    if (magnitude > 1000) points.forEach(point => { point.x = point.x / magnitude * 100; point.y = point.y / magnitude * 100; });
    const minX = Math.min(...points.map(point => point.x)), maxX = Math.max(...points.map(point => point.x)), minY = Math.min(...points.map(point => point.y)), maxY = Math.max(...points.map(point => point.y));
    // Each pair needs separation in just one direction; leave room for node labels and numbers.
    let requiredScale = 90;
    points.forEach((point, i) => points.slice(i + 1).forEach(other => {
      const dx = Math.abs(point.x - other.x), dy = Math.abs(point.y - other.y);
      if (dx < 1e-6 && dy < 1e-6) return;
      requiredScale = Math.max(requiredScale, Math.min(dx > 1e-6 ? 86 / dx : Infinity, dy > 1e-6 ? 74 / dy : Infinity));
    }));
    const scale = Math.min(180, Math.ceil(requiredScale));
    const drawingWidth = Math.max(960, (maxX - minX) * scale + 130), drawingHeight = Math.max(420, (maxY - minY) * scale + 120);
    const placed = [];
    nodes = atoms.map((atom, index) => {
      const origin = { x: drawingWidth / 2 + (points[index].x - (minX + maxX) / 2) * scale, y: drawingHeight / 2 + (points[index].y - (minY + maxY) / 2) * scale };
      const free = point => point.x >= 42 && point.y >= 46 && !placed.some(other => Math.abs(point.x - other.x) < 86 && Math.abs(point.y - other.y) < 74);
      let position = origin;
      // Exact coordinate duplicates occur in some local branches. Move only colliding nodes.
      for (let radius = 1; !free(position); radius++) {
        const candidates = [];
        for (let dx = -radius; dx <= radius; dx++) for (let dy = -radius; dy <= radius; dy++) {
          if (Math.max(Math.abs(dx), Math.abs(dy)) === radius) candidates.push({ x: origin.x + dx * 86, y: origin.y + dy * 74, distance: (dx * 86) ** 2 + (dy * 74) ** 2 });
        }
        candidates.sort((a, b) => a.distance - b.distance || a.y - b.y || a.x - b.x);
        position = candidates.find(free) || origin;
      }
      placed.push(position);
      return { id: index + 1, element: atom.element, hydrogens: atom.hydrogens, charge: atom.charge || 0, label: groupLabel(atom), x: position.x, y: position.y };
    });
    bonds = (compound.bonds || []).map(bond => ({ from: bond.from + 1, to: bond.to + 1, order: bond.order }));
    backbone = (compound.backbone || []).map(index => index + 1); nextId = nodes.length + 1; selected = null; linkStart = null;
    $('[data-ce-formula]').value = compound.condensed || '';
    render(); notify('Struttura caricata. Puoi modificare gruppi, posizioni e legami.');
  }
  function reset() { nodes = []; bonds = []; backbone = []; selected = null; linkStart = null; nextId = 1; $('[data-ce-formula]').value = ''; render(); notify('Il foglio è vuoto. Aggiungi il primo gruppo.'); }
  function click(event) {
    const node = event.target.closest('[data-ce-node]');
    if (node) { useNode(Number(node.dataset.ceNode)); return; }
    const bond = event.target.closest('[data-ce-bond]');
    if (bond) {
      if (mode === 'delete') { const [from, to] = bond.dataset.ceBond.split('-').map(Number); bonds = bonds.filter(item => !(item.from === from && item.to === to)); backbone = []; render(); notify('Legame eliminato.'); }
      else notify('Per eliminare questo legame scegli lo strumento Elimina.');
      return;
    }
    const modeButton = event.target.closest('[data-ce-mode]');
    if (modeButton) {
      mode = modeButton.dataset.ceMode; linkStart = null;
      container.querySelectorAll('[data-ce-mode]').forEach(button => { const active = button.dataset.ceMode === mode; button.classList.toggle('active', active); button.setAttribute('aria-pressed', active); });
      render(); notify({ add: 'Scegli il gruppo e tocca il foglio per aggiungerlo.', link: 'Tocca il primo gruppo, poi il secondo.', move: 'Seleziona il gruppo da spostare, poi tocca la nuova posizione.', delete: 'Tocca il gruppo o il legame da eliminare.' }[mode]); return;
    }
    if (event.target.closest('[data-ce-add]')) { const column = nodes.length % 9, row = Math.floor(nodes.length / 9); add(76 + column * 100, 90 + row * 90); return; }
    if (event.target.closest('[data-ce-reset]')) { reset(); return; }
    if (event.target.closest('[data-ce-connect]')) { connect(Number($('[data-ce-from]').value), Number($('[data-ce-to]').value)); return; }
    if (event.target.closest('[data-ce-delete]')) { const id = Number($('[data-ce-from]').value); if (id) remove(id); else notify('Scegli prima il gruppo da eliminare.', true); return; }
    if (event.target.closest('[data-ce-change]')) {
      const node = nodes.find(item => item.id === Number($('[data-ce-from]').value));
      if (!node) { notify('Scegli prima il gruppo da modificare.', true); return; }
      Object.assign(node, groupAtom()); render(); notify(`Il gruppo ${node.id} ora è ${node.label}.`); return;
    }
    if (event.target.closest('[data-ce-read]')) {
      try { const graph = parseCondensedFormula($('[data-ce-formula]').value); const compound = compoundFromGraph(graph); loadCompound(compound); notify('Formula letta e verificata. Ora puoi modificarla oppure analizzare il disegno.'); }
      catch (error) { notify(error.message || 'La formula non è valida.', true); }
      return;
    }
    if (event.target.closest('[data-ce-analyze]')) {
      try { const compound = compoundFromGraph(drawingGraph()); notify('Struttura verificata: le formule sono state aggiornate.'); onAnalyze(compound); }
      catch (error) { notify(error.message || 'Controlla gruppi e legami del disegno.', true); }
      return;
    }
    if (event.target.closest('[data-ce-canvas]')) {
      const box = canvas.getBoundingClientRect(), x = (event.clientX - box.left) * width / box.width, y = (event.clientY - box.top) * height / box.height;
      if (mode === 'add') add(x, y);
      else if (mode === 'move' && selected !== null) { const node = nodes.find(item => item.id === selected); node.x = Math.max(42, x); node.y = Math.max(46, y); render(); notify('Gruppo spostato.'); }
    }
  }
  function keydown(event) { if (['Enter', ' '].includes(event.key) && event.target.matches('[data-ce-node], [data-ce-bond]')) { event.preventDefault(); click(event); } }
  function change(event) { if (event.target.matches('[data-ce-from]')) { selected = Number(event.target.value) || null; render(); } }
  container.addEventListener('click', click); container.addEventListener('keydown', keydown); container.addEventListener('change', change);
  render();
  return { loadCompound, reset, destroy() { destroyed = true; container.removeEventListener('click', click); container.removeEventListener('keydown', keydown); container.removeEventListener('change', change); container.replaceChildren(); } };
}
