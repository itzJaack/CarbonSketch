// Keep the official Sites workflow unchanged, with ordinary subprocess pipes
// and the installed Git Bash on its session-local PATH for archive packaging.
// Credentials are received through hidden stdin and never written to disk.
import { spawn } from 'node:child_process';
import { existsSync } from 'node:fs';
import path from 'node:path';
const [helper,projectId] = process.argv.slice(2);
if (!helper || !projectId) throw new Error('Provide the official Sites helper path and Site ID.');
const terminal = process.stdin.isTTY;
const input = await new Promise((resolve,reject)=>{
  let data = '';
  const finish = () => {
    process.stdin.removeListener('data',receive);
    process.stdin.removeListener('end',end);
    if (terminal) process.stdin.setRawMode(false);
    process.stdin.pause();
    resolve(data.trim());
  };
  const receive = chunk => {
    data += chunk;
    if (data.includes('\u0003')) { if(terminal)process.stdin.setRawMode(false);reject(new Error('Cancelled')); }
    else if (data.length > 65536) { if(terminal)process.stdin.setRawMode(false);reject(new Error('Input too large')); }
    else if (/[\r\n]/.test(data)) finish();
  };
  const end = ()=>finish();
  if (terminal) process.stdin.setRawMode(true);
  process.stdin.setEncoding('utf8');
  process.stdin.on('data',receive);
  process.stdin.once('end',end);
  process.stderr.write('Ready for Site workflow JSON on stdin (input is hidden).\n');
  process.stdin.resume();
});
try { JSON.parse(input); } catch { throw new Error('Provide JSON on stdin.'); }
const env = {...process.env};
if (process.platform === 'win32') {
  const gitBin = path.join(process.env.ProgramFiles || 'C:\\Program Files','Git','bin');
  if (existsSync(path.join(gitBin,'bash.exe'))) {
    const pathKey = Object.keys(env).find(key=>key.toLowerCase()==='path') || 'Path';
    env[pathKey] = gitBin + path.delimiter + (env[pathKey] || '');
    // GNU tar otherwise interprets a Windows drive-letter archive as a host.
    env.TAR_OPTIONS = ((env.TAR_OPTIONS || '') + ' --force-local').trim();
  }
}
const child = spawn(process.execPath,[helper,'--project-id',projectId],{stdio:['pipe','pipe','pipe'],env});
child.stdout.pipe(process.stdout);
child.stderr.pipe(process.stderr);
child.once('error',()=>{process.stderr.write('Could not launch the official Sites helper.\n');process.exitCode=1;});
child.stdin.on('error',()=>{process.exitCode=1;});
child.stdin.end(input+'\n');
child.once('close',code=>{process.exitCode=code??1;});
