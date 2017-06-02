import os
import vcf
import shutil
import re

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


def splitVcf(filename, output_folder=None, template = None):
    """ Separates a vcf file into indel and snp sections.
        Parameters
        ----------
            filename: path
                THe vcf file to split.
            output_folder: path
                The folder to save the output files to. If None,
                the files will be saved to the folder the original file
                resides.
            output_template: string; default None
                If provided, will be used as the basename for the output files.

        Returns
        -------
            output: dict<>
                * 'indel': path
                    The indel file
                * 'snp': path
                    The snp file
    """

    folder, basename = os.path.split(filename)
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
            copy2(source, destination)
            split_callset[caller_name] = destination
        else:
            result = splitVcf(source, output_folder, **kwargs)
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
    modify_mutect_script = os.path.join(somaticseq_folder, "modify_Mutect.py")
    
    fixed_callset = dict()
    for caller, source in callset.items():
        if 'output_folder' in kwargs:
            output_folder = kwargs['output_folder']
        else:
            output_folder = os.path.dirname(source)

        if 'patientId' in kwargs:
            basename = "{}.{}.corrected.vcf".format(patient_id, caller)
        else:
            basename = os.path.basename(source)
            basename, ext = os.path.basename(basename)
            basename = "{}.corrected.vcf".format(basename)

        destination = os.path.join(output_folder, basename)
        if 'varscan' in caller:
            #modify VJSD . py −method VarScan2 − i n f i l e input . vcf −o u t f i l e output . vcf
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
                infile  = filename,
                outfile = destination)
        else:
            shutil.copy2(source, destination)
        fixed_callset[caller] = destination

    return fixed_callset

if __name__ == "__main__":
    folder = ""

    callset = Callset()
