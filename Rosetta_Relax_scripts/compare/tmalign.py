'''
(c) 2011-2012 Thomas Holder, MPI for Developmental Biology
'''

from __future__ import print_function

__author__ = 'Thomas Holder'
__version__ = '1.1'
__license__ = 'BSD-2-Clause'

from pymol import cmd, CmdException


def save_pdb_without_ter(filename, selection, **kwargs):
    v = cmd.get_setting_boolean('pdb_use_ter_records')
    if v:
        cmd.set('pdb_use_ter_records', 0)
    cmd.save(filename, selection, **kwargs)
    if v:
        cmd.set('pdb_use_ter_records')


def alignwithanymethod(mobile, target, methods='align super cealign tmalign',
                       async_=1, quiet=1, **kwargs):
    import threading
    import time
    methods = methods.split()
    async_, quiet = int(kwargs.pop('async', async_)), int(quiet)
    mobile_obj = cmd.get_object_list('first (' + mobile + ')')[0]

    def myalign(method):
        newmobile = cmd.get_unused_name(mobile_obj + '_' + method)
        cmd.create(newmobile, mobile_obj)
        start = time.time()
        cmd.do('%s mobile=%s in %s, target=%s' % (method, newmobile, mobile, target))
        if not quiet:
            print('Finished: %s (%.2f sec)' % (method, time.time() - start))

    for method in methods:
        if async_:
            t = threading.Thread(target=myalign, args=(method,))
            t.setDaemon(1)
            t.start()
        else:
            myalign(method)


def tmalign(mobile, target, args='', exe='TMalign', ter=0, transform=1, object=None, quiet=0):
    import subprocess
    import tempfile
    import os
    import re

    ter, quiet = int(ter), int(quiet)

    mobile_filename = tempfile.mktemp('.pdb', 'mobile')
    target_filename = tempfile.mktemp('.pdb', 'target')
    matrix_filename = tempfile.mktemp('.txt', 'matrix')
    mobile_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (mobile)
    target_ca_sele = '(%s) and (not hetatm) and name CA and alt +A' % (target)

    if ter:
        save = cmd.save
    else:
        save = save_pdb_without_ter
    save(mobile_filename, mobile_ca_sele)
    save(target_filename, target_ca_sele)

    exe = cmd.exp_path(exe)
    args = [exe, mobile_filename, target_filename, '-m', matrix_filename] + args.split()

    try:
        process = subprocess.Popen(args, stdout=subprocess.PIPE,
                universal_newlines=True)
        lines = process.stdout.readlines()
    except OSError:
        print('Cannot execute "%s", please provide full path to TMscore or TMalign executable' % (exe))
        raise CmdException
    finally:
        os.remove(mobile_filename)
        os.remove(target_filename)

    # TMalign >= 2012/04/17
    if os.path.exists(matrix_filename):
        lines += open(matrix_filename).readlines()
        os.remove(matrix_filename)

    r = None
    re_score = re.compile(r'TM-score\s*=\s*(\d*\.\d*)')
    rowcount = 0
    matrix = []
    line_it = iter(lines)
    headercheck = False
    alignment = []
    for line in line_it:
        if "TM-score=" in line:
           match = re_score.search(line)
           if match is not None:
               r = float(match.group(1))
               return r


def mmalign(mobile, target, args='', exe='MMalign', ter=0, transform=1, quiet=0):
    return tmalign(mobile, target, args, exe, ter, transform, quiet=quiet)

# pymol commands
cmd.extend('alignwithanymethod', alignwithanymethod)
cmd.extend('tmalign', tmalign)

# autocompletion
cmd.auto_arg[0].update({
    'tmalign': cmd.auto_arg[0]['align'],
    'mmalign': cmd.auto_arg[0]['align'],
})
cmd.auto_arg[1].update({
    'tmalign': cmd.auto_arg[1]['align'],
    'mmalign': cmd.auto_arg[1]['align'],
})
