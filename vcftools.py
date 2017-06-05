import os
import vcf
import shutil
import sys

if os.name == 'nt':
    GITHUB_FOLDER = os.path.join(os.getenv('USERPROFILE'), 'Documents', 'Github')
else:
    GITHUB_FOLDER = os.path.join(os.getenv('HOME'), 'Documents', 'Github')
sys.path.append(GITHUB_FOLDER)

import pytools.systemtools as systemtools

def copyVcf(source, destination):
    with open(source, 'r') as input_file:
        reader = vcf.Reader(input_file)
        if 'Varscan' in source:
            reader.formats['DP4'] = reader.formats['DP4']._replace(num=4)
            reader.formats['DP4'] = reader.formats['DP4']._replace(type='Integer')
        with open(destination, 'w') as file1:
            writer = vcf.Writer(file1, reader)
            for record in reader:
                filterOut = '/' in str(record.ALT[0]) or '/' in record.REF
                if not filterOut:
                    try:
                        writer.write_record(record)
                    except ValueError:
                        print(record)
    return destination


def splitVcf(filename, output_folder = None):
    """ Separates a vcf file into indel and snp sections.
        Parameters
        ----------
            filename: string [PATH]
                THe vcf file to split.
            output_folder: path
                The folder to save the output files to. If None,
                the files will be saved to the folder the original file
                resides.

        Returns
        -------
            output: dict<>
                * 'indel': path
                    The indel file
                * 'snp': path
                    The snp file
    """

    path, basename = os.path.split(filename)
    if output_folder is None: output_folder = path
    basename = os.path.basename(basename)
    snp_filename = os.path.join(output_folder, basename + ".snp.vcf")
    indel_filename=os.path.join(output_folder, basename + ".indel.vcf")
    
    with open(filename) as vcf_file:
        reader          = vcf.Reader(vcf_file)

        snp_writer      = vcf.Writer(open(snp_filename,   'w'), reader)
        indel_writer    = vcf.Writer(open(indel_filename, 'w'), reader)
        for record in reader:
            if record.is_snp:
                snp_writer.write_record(record)
            elif record.is_indel:
                indel_writer.write_record(record)
            else:
                print("vcftools.splitVcf: ", record)

    result = {
        'snp': snp_filename,
        'indel': indel_filename

    }
    return result
def splitCallset(callset, output_folder, **kwargs):
    split_callset = dict()
    for caller_name, source in callset.items():
        if 'indel' in caller_name or 'snp' in caller_name:
            destination = os.path.join(output_folder, os.path.basename(source))
            shutil.copy2(source, destination)
            split_callset[caller_name] = destination
        else:
            result = splitVcf(source, output_folder)
            split_callset[caller_name + '-indel'] = result['indel']
            split_callset[caller_name + '-snp']   = result['snp']

    return split_callset
def fixCallerOutputs(callset, somaticseq_folder, **kwargs):
    """
        Required Parameters
        -------------------
            variants: dict<str:str> [dict<caller_name: caller_output>]
                The callset to fix
            somaticseq_folder: str [path]
                Th folder containing the somaticseq program.
        Optional Parameters
        -------------------
            sample: dict<>
            patientId: str

    """


    modify_vjsd_script   = os.path.join(somaticseq_folder, "modify_VJSD.py")
    
    fixed_callset = dict()
    for caller, source in callset.items():
        if 'output_folder' in kwargs:
            output_folder = kwargs['output_folder']
        else:
            output_folder = os.path.dirname(source)

        if 'patientId' in kwargs:
            basename = "{}.{}.corrected.vcf".format(kwargs['patientId'], caller)
        else:
            basename = os.path.basename(source)
            basename, ext = os.path.splitext(basename)
            basename = "{}.corrected.vcf".format(basename)

        destination = os.path.join(output_folder, basename)
        if 'varscan' in caller:
            command = """python3 {program} -method Varscan2 -infile {infile} -outfile {outfile}"""
        elif 'somaticsniper' in caller:
            command = """python3 {program} -method SomaticSniper -infile {infile} -outfile {outfile}"""
        elif 'muse' in caller:
            command = """python3 {program} -method MuSE -infile {infile} -outfile {outfile}"""
        else:
            command = None

        if command:
            command = command.format(
                program = modify_vjsd_script,
                infile  = source,
                outfile = destination)
            systemtools.Terminal(command)
        else:
            shutil.copy2(source, destination)
        fixed_callset[caller] = destination

    return fixed_callset

if __name__ == "__main__":
    pass
