// Windows ConPTY cannot be inherited by some child processes. Keep the official
// Sites workflow unchanged, but give it ordinary pipes for its subprocesses.
// Credentials are received through hidden stdin and never written to disk.
import { spawn } from 'node:child_process';
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
const child = spawn(process.execPath,[helper,'--project-id',projectId],{stdio:['pipe','pipe','pipe'],env:process.env});
child.stdout.pipe(process.stdout);
child.stderr.pipe(process.stderr);
child.once('error',()=>{process.stderr.write('Could not launch the official Sites helper.\n');process.exitCode=1;});
child.stdin.end(input+'\n');
child.once('close',code=>{process.exitCode=code??1;});
