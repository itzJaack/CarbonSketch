import { renderCondensed } from './renderers.js';

const FONT_SIZE = 22;
const SUB_SIZE = 14;
const INK = 'var(--diagram-ink, #253a32)';
const PAPER = 'var(--diagram-paper, #ffffff)';
const ACCENT = 'var(--accent, #338365)';
const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const round = value => Number(value.toFixed(3));
const dimension = (value, fallback) => Number.isFinite(Number(value)) && Number(value) > 0 ? Number(value) : fallback;
const token = (text, script = '') => ({ text: String(text), script });
const plain = tokens => tokens.map(part => part.text).join('');
// Monospaced formulas make the extent calculation independent of differing letter widths.
const textWidth = tokens => tokens.reduce((sum, part) => sum + [...part.text].length * (part.script ? SUB_SIZE : FONT_SIZE) * .6, 0);
const markup = tokens => tokens.map(part => part.script
  ? `<tspan font-size="${SUB_SIZE}" baseline-shift="${part.script}">${escape(part.text)}</tspan>`
  : escape(part.text).replace(/-/g, '–')).join('');
const edgeToken = order => token(order === 2 ? '=' : order === 3 ? '≡' : '-');

/** Only use the chemistry engine's identified parent chain in a connected acyclic graph. */
function identifiedTree(compound) {
  const atoms = compound.atoms || [], bonds = compound.bonds || [], chain = compound.backbone;
  if (!Array.isArray(chain) || !chain.length || new Set(chain).size !== chain.length) return null;
  const heavy = atoms.flatMap((atom, index) => atom.element !== 'H' ? [index] : []);
  if (chain.some(index => !Number.isInteger(index) || !atoms[index] || atoms[index].element === 'H')) return null;
  const graph = atoms.map(() => []), edges = new Set();
  for (const bond of bonds) {
    if (!Number.isInteger(bond.from) || !Number.isInteger(bond.to) || !atoms[bond.from] || !atoms[bond.to] || bond.from === bond.to) return null;
    if (atoms[bond.from].element === 'H' || atoms[bond.to].element === 'H') continue;
    const order = Number(bond.order ?? 1), key = [bond.from, bond.to].sort((a, b) => a - b).join('-');
    if (![1, 2, 3].includes(order) || edges.has(key)) return null;
    edges.add(key);
    graph[bond.from].push({ index: bond.to, order });
    graph[bond.to].push({ index: bond.from, order });
  }
  if (edges.size !== heavy.length - 1) return null;
  const visited = new Set(), pending = [chain[0]];
  while (pending.length) {
    const index = pending.pop();
    if (visited.has(index)) continue;
    visited.add(index);
    graph[index].forEach(next => pending.push(next.index));
  }
  if (visited.size !== heavy.length || chain.some((index, i) => i && !graph[chain[i - 1]].some(next => next.index === index))) return null;
  return { atoms, graph, chain };
}

function atomTokens(atom, index, graph) {
  const given = Number(atom.hydrogens);
  const valence = { C: 4, N: Number(atom.charge) === 1 ? 4 : 3, O: 2, F: 1, Cl: 1, Br: 1, I: 1, S: 2, P: 3 };
  const hydrogens = Number.isFinite(given) && given >= 0 ? given : Math.max(0, (valence[atom.element] || 0) - graph[index].reduce((sum, bond) => sum + bond.order, 0));
  const result = [token(atom.element)];
  if (hydrogens) result.push(token('H'), ...(hydrogens > 1 ? [token(hydrogens, 'sub')] : []));
  const charge = Number(atom.charge || 0);
  if (Number.isFinite(charge) && charge) result.push(token(`${Math.abs(charge) === 1 ? '' : Math.abs(charge)}${charge > 0 ? '+' : '−'}`, 'super'));
  return result;
}

function branchFormula(index, parent, tree, mainSet) {
  const head = atomTokens(tree.atoms[index], index, tree.graph);
  const children = tree.graph[index].filter(next => next.index !== parent && !mainSet.has(next.index));
  const fragments = children.map(next => ({ ...next, fragment: branchFormula(next.index, index, tree, mainSet) }));
  const indices = [index, ...fragments.flatMap(child => child.fragment.indices)];
  if (!fragments.length) return { tokens: head, indices };
  if (fragments.length === 1) return { tokens: [...head, edgeToken(fragments[0].order), ...fragments[0].fragment.tokens], indices };
  // Compact identical rooted subtrees, preserving the bond order of each attachment.
  const groups = new Map();
  for (const child of fragments) {
    const key = `${child.order}:${plain(child.fragment.tokens)}`;
    if (!groups.has(key)) groups.set(key, { child, count: 0 });
    groups.get(key).count++;
  }
  const result = [...head];
  for (const { child, count } of groups.values()) {
    result.push(token('('), ...(child.order === 1 ? [] : [edgeToken(child.order)]), ...child.fragment.tokens, token(')'));
    if (count > 1) result.push(token(count, 'sub'));
  }
  return { tokens: result, indices };
}

function bondLines(from, to, order, attrs = '') {
  const length = Math.hypot(to.x - from.x, to.y - from.y) || 1;
  const nx = -(to.y - from.y) / length, ny = (to.x - from.x) / length;
  const offsets = order === 3 ? [-4.5, 0, 4.5] : order === 2 ? [-2.6, 2.6] : [0];
  return offsets.map(offset => `<line class="bond" data-order="${order}" ${attrs} x1="${round(from.x + nx * offset)}" y1="${round(from.y + ny * offset)}" x2="${round(to.x + nx * offset)}" y2="${round(to.y + ny * offset)}"/>`).join('');
}

function graphFallback(compound, options) {
  const label = escape(`Condensata 2D: formula semisviluppata del grafo, anelli o catena principale non identificata: ${compound.name || 'composto'}`);
  return renderCondensed(compound, options)
    .replace('class="molecule-diagram molecule-condensed"', 'class="molecule-diagram molecule-condensed molecule-condensed2d molecule-condensed2d-fallback" data-fallback="true" data-layout="graph"')
    .replace(/aria-label="[^"]*"/, `aria-label="${label}"`)
    .replace('<desc>', '<desc>Condensata 2D: per questa struttura è usata la formula semisviluppata del grafo; gli anelli e tutti i loro legami restano visibili. ');
}

/** Horizontal condensed formula with graph-derived substituents on their parent carbon. */
export function renderCondensed2D(compound, { width = 1000, height = 380, numbering = false } = {}) {
  width = dimension(width, 1000); height = dimension(height, 380);
  const tree = identifiedTree(compound);
  if (!tree) return graphFallback(compound, { width, height });
  const mainSet = new Set(tree.chain);
  const columns = tree.chain.map((index, position) => {
    const tokens = atomTokens(tree.atoms[index], index, tree.graph), half = textWidth(tokens) / 2;
    const branches = tree.graph[index].filter(next => !mainSet.has(next.index)).map((next, slot) => {
      const fragment = branchFormula(next.index, index, tree, mainSet);
      const rootHalf = textWidth(atomTokens(tree.atoms[next.index], next.index, tree.graph)) / 2;
      const branchWidth = textWidth(fragment.tokens), layer = Math.floor(slot / 2);
      return { ...next, ...fragment, label: plain(fragment.tokens), rootHalf, width: branchWidth, side: slot % 2 ? 1 : -1, dy: (slot % 2 ? 1 : -1) * 88 * (layer + 1), dx: layer ? -(branchWidth + 42) * layer : 0 };
    });
    return { index, position, tokens, label: plain(tokens), half, branches, left: Math.max(half, ...branches.map(branch => branch.rootHalf - branch.dx)), right: Math.max(half, ...branches.map(branch => branch.dx + branch.width - branch.rootHalf)) };
  });
  const margin = 42;
  columns.forEach((column, i) => { column.x = i ? columns[i - 1].x + Math.max(78, columns[i - 1].right + column.left + 28) : margin + column.left; });
  const minimumWidth = columns.at(-1).x + columns.at(-1).right + margin;
  const svgWidth = Math.ceil(Math.max(width, minimumWidth));
  const centerOffset = (svgWidth - minimumWidth) / 2;
  columns.forEach(column => { column.x += centerOffset; });
  const branches = columns.flatMap(column => column.branches);
  const above = Math.max(0, ...branches.filter(branch => branch.side < 0).map(branch => -branch.dy));
  const below = Math.max(0, ...branches.filter(branch => branch.side > 0).map(branch => branch.dy));
  const svgHeight = Math.ceil(Math.max(height, above + below + (numbering && !below ? 72 : 64)));
  const mainY = svgHeight / 2 + (above - below) / 2;
  let body = '<g class="condensed2d-main-bonds">';
  columns.slice(0, -1).forEach((column, i) => {
    const next = columns[i + 1], order = tree.graph[column.index].find(bond => bond.index === next.index).order;
    body += bondLines({ x: column.x + column.half + 8, y: mainY }, { x: next.x - next.half - 8, y: mainY }, order, `data-bond="${column.index}-${next.index}"`);
  });
  body += '</g><g class="condensed2d-main-chain">';
  for (const column of columns) {
    body += `<text class="atom-label main-chain-label" data-atom="${column.index}" data-locant="${column.position + 1}" data-label="${escape(column.label)}" x="${round(column.x - column.half)}" y="${round(mainY)}">${markup(column.tokens)}</text>`;
  }
  body += '</g><g class="condensed2d-branches">';
  for (const column of columns) for (const branch of column.branches) {
    const rootX = column.x + branch.dx, rootY = mainY + branch.dy;
    body += `<g class="condensed2d-branch" data-parent="${column.index}" data-parent-locant="${column.position + 1}" data-root="${branch.index}" data-atoms="${branch.indices.join(',')}" data-side="${branch.side < 0 ? 'above' : 'below'}" data-label="${escape(branch.label)}">`;
    body += bondLines({ x: column.x, y: mainY + branch.side * 19 }, { x: rootX, y: rootY - branch.side * 19 }, branch.order, `data-bond="${column.index}-${branch.index}"`);
    body += `<text class="atom-label branch-label" data-atom="${branch.index}" data-label="${escape(branch.label)}" x="${round(rootX - branch.rootHalf)}" y="${round(rootY)}">${markup(branch.tokens)}</text></g>`;
  }
  body += '</g>';
  if (numbering) {
    body += '<g class="backbone-numbering">';
    for (const column of columns) {
      const x = column.x - (column.branches.some(branch => branch.side > 0) ? column.half + 10 : 0);
      body += `<text class="atom-number" data-atom="${column.index}" data-locant="${column.position + 1}" x="${round(x)}" y="${round(mainY + 30)}">${column.position + 1}</text>`;
    }
    body += '</g>';
  }
  const description = `Condensata 2D: ogni gruppo della catena principale è esplicito, con gli idrogeni del suo atomo. Catena principale: ${columns.map(column => column.label).join('-')}. Sostituenti: ${columns.flatMap(column => column.branches.map(branch => `${column.position + 1}: ${branch.label}`)).join('; ') || 'nessuno'}. I legami multipli sono conservati. La catena principale segue quella identificata dal motore chimico.`;
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${svgWidth} ${svgHeight}" width="${svgWidth}" height="${svgHeight}" class="molecule-diagram molecule-condensed2d" data-layout="horizontal" data-min-width="${Math.ceil(minimumWidth)}" data-main-chain="${tree.chain.join(',')}" role="img" aria-label="${escape(`Formula condensata 2D: ${compound.name || 'composto'}`)}"><title>${escape(compound.name || 'Composto')}</title><desc>${escape(description)}</desc><style>.molecule-condensed2d text{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;dominant-baseline:central;fill:${INK};font-size:${FONT_SIZE}px;font-weight:500}.molecule-condensed2d .atom-label{paint-order:stroke;stroke:${PAPER};stroke-width:4;stroke-linejoin:round}.molecule-condensed2d .bond{stroke:${INK};stroke-width:2.3;stroke-linecap:round}.molecule-condensed2d .atom-number{font-family:Inter,ui-sans-serif,system-ui,sans-serif;fill:${ACCENT};font-size:12px;text-anchor:middle;font-weight:600;paint-order:stroke;stroke:${PAPER};stroke-width:3}</style>${body}</svg>`;
}
