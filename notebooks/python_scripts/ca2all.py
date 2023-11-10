# source: https://bitbucket.org/lcbio/ca2all/src/master/ca2all.py
import os
import glob
import re
import sys
import argparse

from tempfile import mkstemp
from os.path import basename

from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

_PIR_TEMPLATE = '\n'.join(
    ['>P1;%s', 'sequence:::::::::', '%s', '*', '',
        '>P1;model_ca', 'structure:%s:FIRST:@:END:@::::', '*']
)

filename = '/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/data/interim/foldseek_results/pdb_output/IL17B_ENSG00000127743/AF-Q9UHF5-F1-model_v4/AF-Q9UHF5-F1-model_v4_foldseek_Alphafold_UniProt50AF-Q9UHF5-F1-model_v4.pdb_AF-A0A2G6HJI5-F1-model_v4.pdb'
output = "/central/groups/MazmanianLab/joeB/Microbiota-Immunomodulation/Microbiota-Immunomodulation/TEST_ca2all.pdb"
iterations = 3
verbose = True


def ca2all(filename, output=None, iterations=1, verbose=False):

    old_stdout = sys.stdout
    if verbose:
        sys.stdout = sys.stderr
    else:
        sys.stdout = open(os.devnull, 'w')

    pdb = mkstemp(prefix='.', suffix='.pdb', dir='.', text=True)[1]
    prefix = basename(pdb).rsplit('.', 1)[0]
    aa_names = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    aa_names = {v: k for k, v in aa_names.items()}
    atoms = []
    # pattern = re.compile('ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*')
    # pattern = re.compile('ATOM.{10}CA .([A-Z]{3}) ([A-Z]) (.{5}).{27}(.{12}).*')
    pattern = re.compile('ATOM.{10}CA .([A-Z]{3}) ([A-Z])')
    try:
        with open(filename, 'r') as f, open(pdb, 'w') as tmp:
            for line in f:
                if line.startswith('ENDMDL'):
                    break
                else:
                    match = re.match(pattern, line)
                    if match:
                        atoms.append(match.groups())
                        tmp.write(line)
        if not len(atoms):
            raise Exception('File %s contains no CA atoms' % filename)
        chains = [atoms[0][1]]
        seq = ''
        rr = int(atoms[0][2]) - 1
        for a in atoms:
            s, c, r = a[:3]
            if int(r) != int(rr) + 1:
                seq += '/'
            rr = r
            seq += aa_names[s]
            if c not in chains:
                chains += c
        pir = prefix + '.pir'
        with open(pir, 'w') as f:
            f.write(_PIR_TEMPLATE % (prefix, seq, pdb))
##### ---MODELLER SECTION---#####
        env = environ()
        env.io.atom_files_directory = ['.']
        env.libs.topology.read(file='$(LIB)/top_allh.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')

        class MyModel(automodel):
            def special_patches(self, aln):
                self.rename_segments(segment_ids=chains)
        mdl = MyModel(
            env,
            alnfile=pir,
            knowns='model_ca',
            sequence=prefix,
            assess_methods=assess.DOPE
        )
        mdl.md_level = refine.slow
        mdl.auto_align(matrix_file=prefix + '.mat')
        mdl.starting_model = 1
        mdl.ending_model = int(iterations)
        mdl.final_malign3d = True
        mdl.make()
        models = [m for m in mdl.outputs if m['failure'] is None]
        cmp_key = 'DOPE score'
        models.sort(lambda x, y: cmp(x[cmp_key], y[cmp_key]))
        final = models[0]['name'].rsplit('.', 1)[0] + '_fit.pdb'
        sys.stdout = old_stdout
        sss = complete_pdb(env, final)
        sss.write(file='OUT_H.pdb', model_format='PDB')

        if output:
            outfile = open(output, 'w')
        else:
            outfile = sys.stdout
        with open(final) as f:
            a = iter(atoms)
            current = ch = r = t = nl = None
            for line in f:
                if line.startswith('ATOM'):
                    res = line[21:27]
                    if not current or current != res:
                        current = res
                        ch, r, t = a.next()[1:]
                    nl = line[:21] + ch + r + line[27:54] + t
                    if len(line) > 66:
                        nl += line[66:]
                    outfile.write(nl)
                elif line.startswith('TER '):
                    outfile.write(line[:22] + nl[22:27] + '\n')
                else:
                    outfile.write(line)
    finally:
        junk = glob.glob(prefix + '*')
        map(os.remove, junk)


# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(
#         prog='ca2all',
#         description="""
# Program rebuilds all atom representation from CA-only models.
# Models are generated by modeller's automodel method and assessed
# by the DOPE score.""",
#         formatter_class=argparse.RawTextHelpFormatter,
#         epilog=''
#     )
#     parser.add_argument(
#         '-i', '--input-pdb',
#         help='input pdb file with structure in Ca representation',
#         metavar='INPUT',
#         dest='inp',
#         required=True
#     )
#     parser.add_argument(
#         '-o', '--output-pdb',
#         help='save output pdb [default stdout]',
#         default=None,
#         metavar='OUTPUT',
#         dest='out'
#     )
#     parser.add_argument(
#         '-n', '--number-of-iterations',
#         help='number of models generated by modeller',
#         metavar='ITER',
#         type=int,
#         default=1,
#         dest='iter'
#     )
#     parser.add_argument(
#         '-v, --verbose',
#         help='print modeller output to stderr',
#         action='store_true',
#         dest='verbosity'
#     )
#     args = parser.parse_args()
#     ca2all(args.inp, args.out, args.iter, args.verbosity)
