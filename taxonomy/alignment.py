#!/usr/bin/env python
import json, sys, os, string, argparse, subprocess, re
from string import Template
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, SeqRecord


# idea:
# pull the rf line out of the ref sto.
# get the mask in terms of the consensus columns of the ref alignment
# align with hmmalign --mapali
# make the consensus columns alignment
# trim it using the mask


class Alignment(object):
    """
    A class to provide alignment-related tools for reference packages.
    """

    def __init__(self, reference_package, out_prefix, debug=False, verbose=False):
        """
        Constructor - sets up a number of properties when instantiated.
        """

        self.reference_package = reference_package        
        # Determine the name of the reference package, excluding any other 
        # elements in its path.
        self.reference_package_name = list(os.path.split(reference_package)).pop()
        self.reference_package_name_prefix = os.path.splitext(self.reference_package_name)[0]

        # Keep track of whether or not an out_prefix was specifed.
        self.out_prefix_arg = out_prefix
        self.out_prefix = out_prefix
        # Set default prefix if unspecified.
        if out_prefix is None:
            # Sequence file names will be appended within for loops.
            out_prefix = os.path.join(reference_package, reference_package_name_prefix) + '.'

        # Read in CONTENTS.json and retrieve settings we will be working with.
        json_file = os.path.join(reference_package, 'CONTENTS.json')
        json_contents = self._read_contents_json(json_file)
        aln_sto, aln_fasta, profile = [reference_package] * 3
        self.aln_sto = os.path.join(aln_sto, json_contents['files']['aln_sto'])
        self.aln_fasta = os.path.join(aln_fasta, json_contents['files']['aln_fasta'])
        self.profile = os.path.join(profile, json_contents['files']['profile'])
        if 'mask' in json_contents['files']:
            self.mask_file = os.path.join(self.reference_package, json_contents['files']['mask'])

        # read in the consensus RF line
	self.consensus_rf = self._consensus_rf_of_sto(self.aln_sto)
        self.debug = debug
        self.verbose = verbose

            
    # Public methods

    def hmmer_search(self, search_options, sequence_file):
        """
        Recruit fragments using hmmsearch.  Works with a single sequence file, further 
        work would be required if it is to be expanded to work with multiple sequence files.
        """
        # hmmsearch must be in PATH for this to work.
        hmmsearch_output_file = self.out_prefix + '.search_out.sto'
        hmmsearch_command = 'hmmsearch --notextw --noali -A ' + hmmsearch_output_file + \
                            ' ' + self.profile + ' ' + sequence_file


        child = subprocess.Popen(hmmsearch_command,
                                 stdin=None,
                                 stdout=None,
                                 stderr=None,
                                 shell=(sys.platform!="win32"))
        return_code = child.wait()

        # If return code was not 1, hmmsearch completed without errors.
        if not return_code:
            return hmmsearch_output_file
        else:
            raise Exception, "hmmsearch command failed: \n" + hmmsearch_command


    def hmmer_align(self, sequence_files, squeeze=False, mask=False, 
                    frag=False, ref=False, separate_steps=False, sequence_file_format='fasta'):
        """
        Create an alignment with hmmalign. Then, separate out reference sequences 
        from the fragments into two separate files. If separate_steps is True, 
        separate files with squeeze and mask output will be written.  Note that 
        mask=True forces a squeeze action.
        """
        # Check to make sure aln_sto has no all-gap columns.  An exception 
        # is thrown if one or more all-gap columns are found.
        self._validate_stockholm_alignment(self.aln_sto)
       
        # Get first sequence length from an alignment.
        aln_sto_length = self._get_sequence_length(self.aln_sto, 'stockholm')
        aln_fasta_length = self._get_sequence_length(self.aln_fasta, 'fasta')

        # hmmalign must be in PATH for this to work.
        hmmer_template = Template('hmmalign -o $tmp_file' + ' --mapali ' + \
                                  self.aln_sto + ' ' + self.profile + ' $sequence_file')
        for sequence_file in sequence_files:
            
            # Determine a name for the temporary output file.
            tmp_file = sequence_file + '.' + str(os.getpid()) + '.sto'
            _,sequence_file_name = os.path.split(sequence_file)
            sequence_file_name_prefix = string.join(list(os.path.splitext(sequence_file_name))[0:-1])

            hmmalign_command = hmmer_template.substitute(sequence_file=sequence_file,
                                                         aln_sto=self.aln_sto,
                                                         profile=self.profile,
                                                         tmp_file=tmp_file,
                                                        )

            try:
                child = subprocess.Popen(hmmalign_command,
                                         stdin=None,
                                         stdout=None,
                                         stderr=None,
                                         shell=(sys.platform!="win32"))
                return_code = child.wait()

            
                # If return code was not 1, split off the reference sequences from the fragments.
                if not return_code:
                    # Determine output file names.  Set to default if -o was not specified.
                    out_refs, out_frags = [os.path.join(self.out_prefix + '.ref.fasta'), 
                                           os.path.join(self.out_prefix + '.frag.fasta')]
                    if not self.out_prefix_arg:
                        out_refs = self.out_prefix + sequence_file_name_prefix + ".ref.fasta"
                        out_frags = self.out_prefix + sequence_file_name_prefix + ".frag.fasta"

                    frag_names = self._names(SeqIO.parse(sequence_file, sequence_file_format))
   

                    # We need to write out two files, so we need two iterators.
                    in_frags = SeqIO.parse(tmp_file, "stockholm")
                    in_seqs = SeqIO.parse(tmp_file, "stockholm")
     
                    if squeeze or mask_file:
                        # A separate iterator to determine squeeze gaps is necessary
                        squeeze_seqs = SeqIO.parse(tmp_file, "stockholm")
                        squeeze_seqs = self._id_filter(squeeze_seqs, lambda(idstr): idstr not in frag_names)
                        gaps = self._squeeze_gaps(squeeze_seqs)

                        # Setup squeeze generator for in_seqs and in_frags
                        in_seqs = self._squeezerator(in_seqs, gaps)
                        in_frags = self._squeezerator(in_frags, gaps)

                        # Setup sequence length validation generator for in_seqs, 
                        # comparing with aln_sto.
                        in_seqs = self._sequence_length_check(in_seqs, aln_sto_length)

                        if separate_steps:
                            # Write just frag output after only squeezing is finished.
                            if frag:
                                SeqIO.write(self._id_filter(in_frags, lambda(idstr): idstr in frag_names), 
                                            self.out_prefix + '.squeezed.fasta', "fasta")
                                # It is very important to reset in_frags for later separate_steps use.  
                                in_frags = SeqIO.parse(tmp_file, "stockholm")
                                in_frags = self._squeezerator(in_frags, gaps)

                    if mask:
			mask = _mask_of_file(self.mask_file)

                        # Setup mask generator for in_seqs and in_frags
                        in_seqs = self._maskerator(in_seqs, mask)
                        in_frags = self._maskerator(in_frags, mask)

                        # Setup sequence length validation generator for in_seqs, 
                        # comparing with aln_fasta.
                        in_seqs = self._sequence_length_check(in_seqs, aln_fasta_length)

                        if separate_steps:
                            # Write just frag output after both squeezing and masking.
                            if frag:
                                SeqIO.write(self._id_filter(in_frags, lambda(idstr): idstr in frag_names), 
                                            self.out_prefix + '.masked.fasta', "fasta")

                    # Two separate files are written out, so two separate sets of iterators/generators are used.
                    # If separate steps is True, 
                    if ref and not separate_steps:
                        SeqIO.write(self._id_filter(in_seqs, lambda(idstr): idstr not in frag_names), out_refs, "fasta")
                    if frag and not separate_steps:
                        SeqIO.write(self._id_filter(in_frags, lambda(idstr): idstr in frag_names), out_frags, "fasta")
            except:
                raise
            finally:
                # Always remove the temporary alignment file.
                os.remove(tmp_file)


    # Private methods

    def _read_contents_json(self, json_file):
        with open(json_file, 'r') as contents:
            json_contents = json.load(contents)
        return json_contents 


    def _names(self, records):
        """
        Get the sequence names.
        """
        s = set()
        for record in records:
            s.add(record.id)
        return(s)


    def _id_filter(self, in_seqs, f):
        """
        Generator function to filter out sequences.
        """
        for record in in_seqs:
            if (f(record.id)):
                yield record
        
    
    # Begin squeeze-related functions

    def _squeeze_gaps(self, records):
        """
        Determine which gaps can be squeezed out, building up a template.
        """
        # Need to iterate an additional time to determine which 
        # gaps are shared between all sequences in an alignment.
        gaps = []
        for record in records:
            if len(gaps) == 0:
                gaps_length = len(str(record.seq))
                gaps = [1] * gaps_length
            gaps = map(self._gap_check, gaps, list(str(record.seq)))

        return gaps


    def _squeezerator(self, records, gaps):
        """
        Remove any gaps that are present in the same position across all sequences in an alignment.
        """
        sequence_length = len(gaps)
        for record in records:
            sequence = list(str(record.seq))
            squeezed = []
            position = 0
            while (position < sequence_length):
                if bool(gaps[position]) is False:
                    squeezed.append(sequence[position])
                position += 1
            yield SeqRecord(Seq(''.join(squeezed)), id=record.id,
                            description=record.description)


    def _is_gap(self, character):
        """
        Find out if the current position has a gap.  '.' is 
        used to represent gaps in the stockholm format, but Biopython 
        SeqRecords seem to store these as '-'.
        """
        if character == '-':
            return 1
        else:
            return 0


    def _gap_check(self, gap, character):
        """
        Build up a gaps list that is used on all sequences 
        in an alignment.
        """
        # Skip any characters that have already been found
        if gap == 0:
            return gap
        return int(bool(gap) & bool(self._is_gap(character)))

    # End squeeze-related functions


    # Begin mask-related functions
    def _maskerator(self, records, mask):
        """
        Prune down all sequences to a select list of positions, a mask.
        """
        for record in records:
            sequence = list(str(record.seq))
            yield SeqRecord(Seq(''.join([sequence[i] for i in mask])), 
                            id=record.id, description=record.description)

    def _mask_of_file(self, mask_file):
        """
        Get an integer list from a file.
        """
        # Regex to remove whitespace from mask file.
        whitespace = re.compile(r'\s|\n', re.MULTILINE)
        with open(mask_file, 'r') as handle:
            mask_text = handle.read() 
            mask_text = re.sub(whitespace, '', mask_text)

        # Cast mask positions to integers
        return(map(int, mask_text.split(',')))


    # End mask-related functions


    # Alignment-validation-related functions

    def _validate_stockholm_alignment(self, aln_sto):
        """
        Make sure there are no all-gap columns in aln_sto.
        """
        align_sto = AlignIO.read(aln_sto, 'stockholm')
        gaps = self._squeeze_gaps(align_sto)
        if 1 in gaps:
            raise Exception, 'file ' + aln_sto + ' contains ' + \
                             str(gaps.count(1)) + ' all-gap column(s).'

        return 


    def _get_sequence_length(self, source_file, source_file_type):
        """
        Returns the length of the first sequence in an alignment.  If file is empty 
        or cannot be read in by AlignIO.read, an exception will be thrown.
        """
        alignment = AlignIO.read(source_file, source_file_type)
        return len(alignment.__getitem__(0))


    def _sequence_length_check(self, records, reference_length):
        """
        Generator to validate sequence lengths, as we go through the sequences.
        """
        for record in records:
            if reference_length == len(record):
                yield record
            else:
                raise Exception, 'Sequence "' + record.id + '" has length of ' + str(len(record)) + \
                                 ', expected ' + str(reference_length)

    # End alignment-validation-related functions


    # consensus column related functions

    def _consensus_rf_of_sto(self, sto_aln):
	rf_rex = re.compile("#=GC RF\s+([x.]*)")
	rf_list = []
	with open(sto_aln, 'r') as handle:
	    for line in handle.readlines():
		m = rf_rex.match(line)
		if m:
		    rf_list.append(m.group(1))
            return("".join(rf_list))


    # End consensus column related functions

