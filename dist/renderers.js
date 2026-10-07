import { buildLewisGraph } from './chemistry.js';

const COLORS = { C: '#354b43', H: '#c7d3ce', O: '#dd6875', N: '#688dcc', F: '#72b59b', Cl: '#72b59b', Br: '#b7835f', I: '#ad82bc', S: '#d5b04d', P: '#d8925a' };
const INK = 'var(--diagram-ink, #253a32)';
const MUTED = 'var(--diagram-muted, #71857d)';
const PAPER = 'var(--diagram-paper, #ffffff)';
const ACCENT = 'var(--accent, #338365)';
const TAU = Math.PI * 2;
const finite = value => Number.isFinite(Number(value));
const round = value => Number(value.toFixed(3));
const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

function adjacency(atoms, bonds) {
  const result = atoms.map(() => []);
  bonds.forEach(bond => {
    if (!result[bond.from] || !result[bond.to]) return;
    result[bond.from].push({ index: bond.to, order: bond.order || 1 });
    result[bond.to].push({ index: bond.from, order: bond.order || 1 });
  });
  return result;
}

// This fallback is a two-dimensional molecular diagram, independent of the 3D conformer.
function fallbackLayout(compound) {
  const atoms = compound.atoms || [];
  const neighbors = adjacency(atoms, compound.bonds || []);
  let chain = Array.isArray(compound.backbone) ? compound.backbone.filter(i => atoms[i]) : [];
  if (!chain.length && atoms.length) {
    let best = [];
    function visit(index, path) {
      if (path.length > best.length) best = path;
      if (path.length > 48) return;
      for (const next of neighbors[index]) if (!path.includes(next.index)) visit(next.index, [...path, next.index]);
    }
    // Restrict traversal work for large remote molecules; their PubChem 2D coordinates are normally present.
    if (atoms.length <= 24 && compound.bonds.length < atoms.length + 3) atoms.forEach((_, i) => visit(i, [i]));
    chain = best.length ? best : [0];
  }
  const points = new Array(atoms.length);
  chain.forEach((index, i) => { points[index] = { x: i * 1.35, y: i % 2 ? -0.42 : 0.42 }; });
  function branch(index, parent, angle, depth = 0) {
    if (points[index] || depth > atoms.length) return;
    const origin = points[parent];
    points[index] = { x: origin.x + Math.cos(angle) * 1.45, y: origin.y + Math.sin(angle) * 1.45 };
    const children = neighbors[index].filter(next => !points[next.index]);
    children.forEach((next, j) => branch(next.index, index, angle + (children.length > 1 ? (j - (children.length - 1) / 2) * 1.2 : depth % 2 ? 0.48 : -0.48), depth + 1));
  }
  chain.forEach((index, i) => {
    const children = neighbors[index].filter(next => !points[next.index]);
    const side = i % 2 ? -1 : 1;
    children.forEach((next, j) => branch(next.index, index, side * Math.PI / 2 + (j - (children.length - 1) / 2) * 1.1));
  });
  atoms.forEach((_, i) => { if (!points[i]) points[i] = { x: i * 1.4, y: 0 }; });
  return points;
}

function rawPositions(compound) {
  if ((compound.atoms || []).every(atom => finite(atom.x2d) && finite(atom.y2d))) {
    return compound.atoms.map(atom => ({ x: Number(atom.x2d), y: -Number(atom.y2d) }));
  }
  return fallbackLayout(compound);
}

function fitPositions(points, width, height, padding, maxScale) {
  if (!points.length) return [];
  const xs = points.map(p => p.x), ys = points.map(p => p.y);
  const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
  const spanX = Math.max(maxX - minX, 0.8), spanY = Math.max(maxY - minY, 0.8);
  const scale = Math.min((width - padding * 2) / spanX, (height - padding * 2) / spanY, maxScale);
  const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
  return points.map(point => ({ x: width / 2 + (point.x - cx) * scale, y: height / 2 + (point.y - cy) * scale }));
}

function svgShell(compound, mode, width, height, body, description) {
  const labels = { lewis: 'Struttura di Lewis', topology: 'Formula topologica', condensed: 'Formula semisviluppata' };
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" class="molecule-diagram molecule-${mode}" role="img" aria-label="${escape(`${labels[mode] || 'Struttura'}: ${compound.name || 'composto'}`)}"><title>${escape(compound.name || 'Composto')}</title><desc>${escape(description)}</desc><style>.molecule-diagram{width:100%;height:100%;overflow:visible}.molecule-diagram text{font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}.molecule-diagram .atom-label{paint-order:stroke;stroke:${PAPER};stroke-width:6;stroke-linejoin:round}.molecule-diagram .bond{stroke:${INK};stroke-width:2.6;stroke-linecap:round}.molecule-diagram .atom-number{fill:${ACCENT};font-size:12px;font-weight:600;paint-order:stroke;stroke:${PAPER};stroke-width:4}.molecule-diagram .electron{fill:${INK}}</style>${body}</svg>`;
}

function bondSvg(bond, from, to, fromInset = 0, toInset = 0, gap = 5) {
  const dx = to.x - from.x, dy = to.y - from.y, length = Math.hypot(dx, dy) || 1;
  fromInset = Math.min(fromInset, length * 0.3);
  toInset = Math.min(toInset, length * 0.3);
  const ux = dx / length, uy = dy / length, nx = -uy, ny = ux;
  const order = Number(bond.order || 1);
  const unknown = ![1, 1.5, 2, 3].includes(order);
  const offsets = order === 3 ? [-gap, 0, gap] : order === 2 || order === 1.5 ? [-gap / 2, gap / 2] : [0];
  return offsets.map((offset, i) => {
    const aromatic = order === 1.5 && i === 1;
    const x1 = from.x + ux * fromInset + nx * offset, y1 = from.y + uy * fromInset + ny * offset;
    const x2 = to.x - ux * toInset + nx * offset, y2 = to.y - uy * toInset + ny * offset;
    return `<line class="bond${aromatic ? ' aromatic' : unknown ? ' unknown-bond' : ''}" data-bond="${bond.from}-${bond.to}" data-order="${order}" x1="${round(x1)}" y1="${round(y1)}" x2="${round(x2)}" y2="${round(y2)}"${aromatic || unknown ? ' stroke-dasharray="4 5"' : ''}${unknown ? ' opacity="0.65"' : ''}/>`;
  }).join('');
}

function freeAngles(index, points, neighbors) {
  const angles = (neighbors[index] || []).map(next => Math.atan2(points[next.index].y - points[index].y, points[next.index].x - points[index].x)).map(angle => (angle + TAU) % TAU).sort((a, b) => a - b);
  if (!angles.length) return [{ angle: -Math.PI / 2, size: TAU }];
  return angles.map((angle, i) => {
    const next = i === angles.length - 1 ? angles[0] + TAU : angles[i + 1];
    return { angle: (angle + next) / 2, size: next - angle };
  }).sort((a, b) => b.size - a.size);
}

function chargeLabel(charge) {
  const n = Number(charge || 0);
  return n ? `${Math.abs(n) === 1 ? '' : Math.abs(n)}${n > 0 ? '+' : '−'}` : '';
}

function implicitHydrogens(atom, index, neighbors) {
  if (finite(atom.hydrogens)) return Math.max(0, Number(atom.hydrogens));
  const valences = { C: 4, N: Number(atom.charge) === 1 ? 4 : 3, O: 2, F: 1, Cl: 1, Br: 1, I: 1, S: 2, P: 3 };
  return Math.max(0, (valences[atom.element] || 0) - neighbors[index].reduce((sum, bond) => sum + bond.order, 0));
}

export function renderTopology(compound, { numbering = false, width = 760, height = 380 } = {}) {
  const atoms = compound.atoms || [], bonds = compound.bonds || [];
  const neighbors = adjacency(atoms, bonds);
  const points = fitPositions(rawPositions(compound), width, height, 54, 63);
  const labelled = atoms.map(atom => atom.element !== 'C' || !!atom.charge || atoms.length === 1);
  let body = '<g class="molecular-bonds">';
  bonds.forEach(bond => {
    if (!points[bond.from] || !points[bond.to] || atoms[bond.from].element === 'H' || atoms[bond.to].element === 'H') return;
    body += bondSvg(bond, points[bond.from], points[bond.to], labelled[bond.from] ? 13 : 0, labelled[bond.to] ? 13 : 0);
  });
  body += '</g><g class="molecular-atoms">';
  atoms.forEach((atom, index) => {
    if (!labelled[index] || atom.element === 'H') return;
    const point = points[index], hydrogens = implicitHydrogens(atom, index, neighbors);
    const hydrogenText = hydrogens ? `H${hydrogens > 1 ? `<tspan font-size="12" baseline-shift="sub">${hydrogens}</tspan>` : ''}` : '';
    const label = `${escape(atom.element)}${hydrogenText}${chargeLabel(atom.charge) ? `<tspan font-size="11" baseline-shift="super">${escape(chargeLabel(atom.charge))}</tspan>` : ''}`;
    body += `<text class="atom-label" data-atom="${index}" x="${round(point.x)}" y="${round(point.y)}" fill="${INK}" font-size="20" font-weight="500" text-anchor="middle" dominant-baseline="central">${label}</text>`;
  });
  body += '</g>';
  if (numbering) {
    // Number only an explicitly identified parent carbon chain, never PubChem's arbitrary atom order.
    const chain = compound.backbone?.filter(index => atoms[index]?.element === 'C') || [];
    body += '<g class="backbone-numbering">';
    chain.forEach((index, position) => {
      if (!points[index]) return;
      const free = freeAngles(index, points, neighbors)[0].angle;
      body += `<text class="atom-number" data-atom="${index}" data-locant="${position + 1}" x="${round(points[index].x + Math.cos(free) * 23)}" y="${round(points[index].y + Math.sin(free) * 23)}" text-anchor="middle" dominant-baseline="central">${position + 1}</text>`;
    });
    body += '</g>';
  }
  return svgShell(compound, 'topology', width, height, body, 'Ogni vertice e ogni estremità senza etichetta rappresenta un atomo di carbonio; gli idrogeni legati al carbonio sono impliciti. I legami multipli sono mostrati con linee parallele.');
}

/** Semi-developed formula: H are grouped with their own atoms and every ring bond is retained. */
export function renderCondensed(compound, { width = 760, height = 380 } = {}) {
  const atoms = compound.atoms || [], bonds = compound.bonds || [];
  const neighbors = adjacency(atoms, bonds);
  const points = fitPositions(rawPositions(compound), width, height, 58, 68);
  const hydrogens = atoms.map((atom, index) => implicitHydrogens(atom, index, neighbors));
  const insets = atoms.map((atom, index) => hydrogens[index] ? 21 : atom.element.length > 1 ? 15 : 11);
  let body = '<g class="molecular-bonds">';
  bonds.forEach(bond => {
    if (!points[bond.from] || !points[bond.to] || atoms[bond.from].element === 'H' || atoms[bond.to].element === 'H') return;
    body += bondSvg(bond, points[bond.from], points[bond.to], insets[bond.from], insets[bond.to]);
  });
  body += '</g><g class="molecular-atoms">';
  atoms.forEach((atom, index) => {
    if (atom.element === 'H') return;
    const point = points[index], h = hydrogens[index];
    const label = `${escape(atom.element)}${h ? `H${h > 1 ? `<tspan font-size="12" baseline-shift="sub">${h}</tspan>` : ''}` : ''}${chargeLabel(atom.charge) ? `<tspan font-size="11" baseline-shift="super">${escape(chargeLabel(atom.charge))}</tspan>` : ''}`;
    body += `<text class="atom-label" data-atom="${index}" data-element="${escape(atom.element)}" data-hydrogens="${h}" x="${round(point.x)}" y="${round(point.y)}" fill="${INK}" font-size="20" font-weight="500" text-anchor="middle" dominant-baseline="central">${label}</text>`;
  });
  body += '</g>';
  return svgShell(compound, 'condensed', width, height, body, 'Formula semisviluppata: gli idrogeni sono raggruppati con ciascun atomo cui sono legati. Tutti i legami tra atomi non idrogeno, comprese le chiusure degli anelli, restano espliciti.');
}

export function renderLewis(compound, { width = 760, height = 420 } = {}) {
  const graph = buildLewisGraph(compound);
  const atoms = graph.atoms || [], bonds = graph.bonds || [];
  const points = fitPositions(rawPositions(graph), width, height, 48, 58);
  const neighbors = adjacency(atoms, bonds);
  const minBond = bonds.length ? Math.min(...bonds.map(bond => Math.hypot(points[bond.from].x - points[bond.to].x, points[bond.from].y - points[bond.to].y))) : 34;
  // Long chains keep explicit hydrogens without reversed or hidden short bonds.
  const symbolScale = Math.min(1, Math.max(0.38, minBond / 34));
  let body = '<g class="molecular-bonds">';
  bonds.forEach(bond => {
    if (!points[bond.from] || !points[bond.to]) return;
    body += bondSvg(bond, points[bond.from], points[bond.to], 11 * symbolScale, 11 * symbolScale, 4.8 * symbolScale);
  });
  body += '</g><g class="molecular-atoms">';
  atoms.forEach((atom, index) => {
    const point = points[index], h = atom.element === 'H';
    const fill = h ? MUTED : atom.element === 'C' ? INK : COLORS[atom.element] || INK;
    body += `<text class="atom-label" data-atom="${index}" data-element="${escape(atom.element)}" x="${round(point.x)}" y="${round(point.y)}" fill="${fill}" font-size="${round((h ? 16 : 19) * symbolScale)}" style="stroke-width:${round(6 * symbolScale)}" font-weight="${h ? 400 : 600}" text-anchor="middle" dominant-baseline="central">${escape(atom.element)}${chargeLabel(atom.charge) ? `<tspan font-size="${round(11 * symbolScale)}" baseline-shift="super">${escape(chargeLabel(atom.charge))}</tspan>` : ''}</text>`;
    const count = Number(atom.lonePairs || 0);
    if (!count) return;
    const spaces = freeAngles(index, points, neighbors);
    const widest = spaces[0];
    for (let pair = 0; pair < count; pair++) {
      // Divide the largest free sector when it can hold several pairs; otherwise use distinct sectors.
      const angle = widest.size > Math.PI * 1.25 || spaces.length < count
        ? widest.angle + (pair - (count - 1) / 2) * Math.min(0.68, widest.size / (count + 1))
        : spaces[pair % spaces.length].angle;
      const cx = point.x + Math.cos(angle) * 20 * symbolScale, cy = point.y + Math.sin(angle) * 20 * symbolScale;
      const dx = -Math.sin(angle) * 2.5 * symbolScale, dy = Math.cos(angle) * 2.5 * symbolScale;
      body += `<g class="lone-pair" data-atom="${index}"><circle class="electron" cx="${round(cx + dx)}" cy="${round(cy + dy)}" r="${round(1.55 * symbolScale)}"/><circle class="electron" cx="${round(cx - dx)}" cy="${round(cy - dy)}" r="${round(1.55 * symbolScale)}"/></g>`;
    }
  });
  body += '</g>';
  const notes = graph.notes?.length ? ` ${graph.notes.join(' ')}` : '';
  return svgShell(compound, 'lewis', width, height, body, `Tutti gli atomi di carbonio e idrogeno sono espliciti. Ogni linea rappresenta una coppia di elettroni condivisa e ogni coppia di punti una coppia non condivisa.${notes}`);
}

export function getLewisNotes(compound) {
  return buildLewisGraph(compound).notes || [];
}

function molAtoms(compound) {
  const graph = buildLewisGraph(compound);
  const atoms = graph.atoms.map((atom, index) => ({ serial: index + 1, elem: atom.element, x: Number(atom.x) || 0, y: Number(atom.y) || 0, z: Number(atom.z) || 0, bonds: [], bondOrder: [] }));
  graph.bonds.forEach(bond => {
    if (!atoms[bond.from] || !atoms[bond.to]) return;
    atoms[bond.from].bonds.push(bond.to);
    atoms[bond.from].bondOrder.push(Number(bond.order || 1));
    atoms[bond.to].bonds.push(bond.from);
    atoms[bond.to].bondOrder.push(Number(bond.order || 1));
  });
  return atoms;
}

export function create3DViewer(container, compound, options = {}) {
  if (!container) throw new Error('Il contenitore 3D non è disponibile.');
  let viewer, observer, destroyed = false, initialView;
  let style = options.style || 'ball-stick', hydrogens = options.hydrogens !== false;
  const empty = { available: false, setStyle() {}, setHydrogens() {}, setSpin() {}, reset() {}, resize() {}, destroy() { container.replaceChildren(); }, exportPng() { return null; } };
  try {
    if (!globalThis.$3Dmol) throw new Error('Il motore 3D non è stato caricato.');
    container.replaceChildren();
    viewer = globalThis.$3Dmol.createViewer(container, { backgroundColor: '#f5f8f5', backgroundAlpha: 0, antialias: true });
    viewer.setBackgroundColor('#f5f8f5', 0);
    if (compound.sdf3d) viewer.addModel(compound.sdf3d, 'sdf', { keepH: true });
    else viewer.addModel().addAtoms(molAtoms(compound));

    function applyStyle() {
      if (destroyed) return;
      const palette = { C: '#4f7466', H: '#d5e1db', O: '#e37482', N: '#789fda', F: '#77bf9c', Cl: '#77bf9c', Br: '#b68b66', I: '#b891c5', S: '#dfbc55', P: '#e5a36b' };
      const colorscheme = { prop: 'elem', map: palette };
      const settings = style === 'spacefill'
        ? { sphere: { scale: 0.86, colorscheme } }
        : style === 'stick'
          ? { stick: { radius: 0.19, colorscheme } }
          : { stick: { radius: 0.14, colorscheme }, sphere: { scale: 0.26, colorscheme } };
      viewer.setStyle({}, settings);
      if (!hydrogens) viewer.setStyle({ elem: 'H' }, { sphere: { hidden: true }, stick: { hidden: true }, line: { hidden: true } });
      viewer.render();
    }

    applyStyle();
    viewer.zoomTo();
    viewer.rotate(18, 'x');
    viewer.rotate(-12, 'y');
    viewer.render();
    initialView = viewer.getView().slice();
    const resize = () => {
      if (destroyed || !container.clientWidth || !container.clientHeight) return;
      viewer.resize();
      viewer.render();
    };
    if (typeof ResizeObserver !== 'undefined') {
      observer = new ResizeObserver(resize);
      observer.observe(container);
    }
    const canvas = container.querySelector('canvas');
    if (canvas) {
      canvas.setAttribute('aria-label', `Modello tridimensionale interattivo di ${compound.name || 'composto'}. Trascina per ruotare e usa la rotella o due dita per lo zoom.`);
      canvas.setAttribute('role', 'img');
    }
    return {
      available: true,
      setStyle(value) { style = ['ball-stick', 'spacefill', 'stick'].includes(value) ? value : 'ball-stick'; applyStyle(); },
      setHydrogens(value) { hydrogens = !!value; applyStyle(); },
      setSpin(value) { if (!destroyed) viewer.spin(value ? 'y' : false, 0.35, true); },
      reset() { if (!destroyed) { viewer.setView(initialView.slice()); viewer.render(); } },
      resize,
      exportPng() { if (destroyed) return null; viewer.render(); return viewer.pngURI(); },
      destroy() {
        if (destroyed) return;
        destroyed = true;
        observer?.disconnect();
        viewer.spin(false);
        viewer.removeAllModels();
        viewer.removeAllLabels();
        viewer.removeAllShapes();
        const canvases = Array.from(container.querySelectorAll('canvas'));
        canvases.forEach(item => {
          try { const gl = item.getContext('webgl2') || item.getContext('webgl') || item.getContext('experimental-webgl'); gl?.getExtension('WEBGL_lose_context')?.loseContext(); } catch {}
        });
        container.replaceChildren();
      }
    };
  } catch (error) {
    observer?.disconnect();
    if (viewer) { try { viewer.spin(false); viewer.removeAllModels(); } catch {} }
    container.replaceChildren();
    const message = document.createElement('div');
    message.className = 'viewer-error';
    message.textContent = 'La vista 3D richiede WebGL. Riattiva l’accelerazione grafica del browser per esplorare la molecola.';
    container.appendChild(message);
    return { ...empty, error: error.message };
  }
}
