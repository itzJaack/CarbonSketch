import http from 'node:http';
import { readFile } from 'node:fs/promises';
import path from 'node:path';
const root = path.resolve('dist');
const types = { '.html':'text/html; charset=utf-8', '.css':'text/css; charset=utf-8', '.js':'text/javascript; charset=utf-8', '.svg':'image/svg+xml', '.json':'application/json', '.woff2':'font/woff2' };
const port = Number(process.env.PORT || 4173);
http.createServer(async (req, res) => {
  try {
    const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname);
    const filename = path.resolve(root, '.' + (pathname === '/' ? '/index.html' : pathname));
    if (!filename.startsWith(root + path.sep)) { res.writeHead(403).end(); return; }
    const data = await readFile(filename);
    res.writeHead(200, {'Content-Type':types[path.extname(filename)] || 'application/octet-stream', 'Cache-Control':'no-store'}).end(data);
  } catch { res.writeHead(404).end('File non trovato'); }
}).listen(port, '127.0.0.1', () => console.log(`CarbonSketch: http://127.0.0.1:${port}`));
