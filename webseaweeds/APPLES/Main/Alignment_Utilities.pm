### (c) copyright University of Warwick 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Alignment_Utilities Class ###
# class abstracting alignment algorithms, may expand in future to support WU-BLAST, PatternHunter, or other alignment algorithms
use MooseX::Declare;

class Alignment_Utilities {
    use Genome_DB_Utilities;
    use General_Utilities;
    use File::Temp qw (tempdir tempfile);
    use APPLES_Datatypes qw (SequenceFunctionalType TranscriptChoice);
    use Job_Handler;
    use APPLES_Datatypes qw (Boolean);
    use constant {FALSE => 0,
		  TRUE	=> 1};
    use Cwd;
    use Bio::SearchIO;
    use Bio::SearchIO::FastHitEventBuilder;
    use Data::Dumper;

    my $GU = General_Utilities->new();

    method align_protein_sequences_in_genome_DBs (Genome_Sequence_Database_Parameters $query_database, Genome_Sequence_Database_Parameters $target_database, Alignment_Parameters $alignment_parameters, Int $max_number_of_matches, TranscriptChoice $transcript_choice) {
      ### NONE OF THIS IS BEING CACHED/STORED ###
      if ($alignment_parameters->isa('NCBI_BLAST_Parameters')) {
	my $type = 'protein';   
	# make query FASTA file (temporary)
	my $queryfastafile = $self->private_get_fasta_file_for_proteins($query_database, $transcript_choice);
	$GU->user_info(1, "query fasta file = ".$queryfastafile."\n");
	# make target FASTA file (temporary) and format into a BLAST-database
	my $blast_db_name = $self->private_formatdb_ncbi($target_database, $type, $queryfastafile, $transcript_choice);
	$GU->user_info(1, "blast database name = ".$blast_db_name."\n");
	# call BLAST
	my $blast_result = $self->private_blast_fasta_file_against_blast_database($queryfastafile, $blast_db_name, $alignment_parameters);
	# process results
	unlink $queryfastafile; # uncomment for debugging - examine this file
	my $result = $self->private_get_sequence_to_sequence_matchings_from_m8_file($blast_result,$max_number_of_matches);
	# delete BLAST-result
	unlink $blast_result; # comment this line out for debugging
	return $result; # Alignment_Set object containing hashref
      } else {	    
	die 'method not yet implemented for this type of alignment parameters.'
      }
    } # align_protein_sequences_in_genome_DBs #
    
    method private_get_fasta_file_for_proteins (Genome_Sequence_Database_Parameters $database, TranscriptChoice $transcript_choice) {
	my $tempdir = $GU->get_temp_random_directory(FALSE);
	my $tempfile = $tempdir.'/protein_sequences.fasta';
	my $gdbu = Genome_DB_Utilities->new();
	$gdbu->proteins_to_fasta_file($database, $tempfile, $transcript_choice);
	return $tempfile;
    } # private_get_fasta_file_for_proteins #

    method private_get_fasta_file_for_genomic_dna (Genome_Sequence_Database_Parameters $database) {
	my $tempdir = $GU->get_temp_random_directory(FALSE);
	my $tempfile = $tempdir.'/genomic_dna_sequences.fasta';
	my $gdbu = Genome_DB_Utilities->new();
	$gdbu->genome_to_fasta_file($database, $tempfile);
	return $tempfile;
    } # private_get_fasta_file_for_genomic_dna #

    method private_determine_blast_db_name (Genome_Sequence_Database_Parameters $database, SequenceFunctionalType $type, Str $queryfastafile) {
	my $db_type_specific_prefix;
	if ($database->isa('FASTA_Sequence_Database_Parameters')) {
	    my $md5sum = $GU->md5sum($queryfastafile);
	    $db_type_specific_prefix = 'fastadb_'.$md5sum.'_';
	}
	elsif ($database->isa('Ensembl_Database_Parameters')) {
	    $db_type_specific_prefix = 'ensembldb_';
	}
	elsif ($database->isa('Genbank_Sequence_Database_Parameters')) {
	    $db_type_specific_prefix = 'genbankdb_';
	}
	else {
	    die 'prefix for BLAST-DB name not defined yet for this type of sequence database.';
	}
	my $result = $db_type_specific_prefix.$type.'_'.$database->dbname;
	return $result;
    } # private_determine_blast_db_name #

    method private_formatdb_ncbi (Genome_Sequence_Database_Parameters $database, SequenceFunctionalType $type, Str $queryfastafile, TranscriptChoice $transcript_choice) {
	# turns a FASTA file into a BLAST database
	my $gdbu = Genome_DB_Utilities->new();
	my $BLASTbasename = $self->private_determine_blast_db_name($database, $type, $queryfastafile);
	$GU->user_info(3, "BLAST base name = ".$BLASTbasename."\n");
	my $dir = $ENV{'blast_db'};
	$GU->user_info(3, "blast_db directory: ". $dir."\n");
	my $fulldbpath = $dir;
	$GU->user_info(3, "full blast database path: ". $fulldbpath."\n");
	my $fulldbname = $fulldbpath.$BLASTbasename;
	$GU->user_info(3, "full database name: ". $fulldbname."\n");
	if (-e $fulldbname) {
	    # not employed Job_Handler here as it is not set up to regard files as method results
	    # (it only picks up the possibly altered object and the returned data)
	  $GU->user_info(3,$fulldbname." exists\n");
	    return $fulldbname;
	}
	else {
	    my $cmd;
	    my $targetfastafile;
	    my $filecheck = TRUE;
	    if ( $type eq 'genomic_dna' ) {
		$targetfastafile = $self->private_get_fasta_file_for_genomic_dna($database);
		my $character_count = $gdbu->count_number_of_characters_fasta_file($targetfastafile);
		my $dna_check = $gdbu->check_dna_fasta_file($targetfastafile);
		$cmd = 'formatdb -pF -oT -i '. $targetfastafile;
		if (!$dna_check||($character_count==0)) {	
		    $filecheck = FALSE;
		}	
	    }
	    elsif ( $type eq 'protein' ) {
	      $GU->user_info(3, "making target blast database from protein sequences\n");
	      $GU->user_info(3, "making FASTA file\n");
	      $targetfastafile = $self->private_get_fasta_file_for_proteins($database, $transcript_choice);
	      $GU->user_info(3, "checking FASTA file contents\n");
	      my $character_count = $gdbu->count_number_of_characters_fasta_file($targetfastafile);
	      my $protein_check = $gdbu->roughly_check_protein_fasta_file($targetfastafile);
	      $cmd = '/common/bin/formatdb -pT -oT -i '. $targetfastafile; # needs to take /common/bin setting from APPLES.dat via the load_includes method
	      $GU->user_info(3, $cmd."\n");
	      if (!$protein_check||($character_count==0)) {	
		$filecheck = FALSE;		
	      }
	    }
	    else {
	      die "support for this functional type of sequence not implemented yet.";
	    }
	    if (!$filecheck) {
		unlink $targetfastafile;
		die "sequence does not meet requirements for creation of BLAST-DB.";
	    }
	    $cmd = $cmd." -n ". $BLASTbasename."\n";
	    $GU->user_info( 1, "formatting FASTA file into BLAST database ...\n" );
	    my $prev_dir = getcwd();
	    chdir($fulldbpath);
	    # check for $targetfastafile existence:
	    if (-e $targetfastafile) {
	      $GU->user_info(1, $targetfastafile. " EXISTS!\n");
	    }
	    else {
	      $GU->user_info(1, $targetfastafile. " does not exist (yet)!\n");
	    }
	    $GU->user_info(3, "issuing command: ".$cmd." in ".$fulldbpath."\n");
	    
	    system($cmd) == 0 or die "could not execute this command $cmd $! $?\n"; # this is failing with 'Can't exec "formatdb": No such file or directory at' - so does the input file exist yet? YES!
	    chdir($prev_dir);
	    unlink $targetfastafile; # comment this line for debugging
	}
	return($fulldbname);
    } # private_formatdb_ncbi #

    method private_blast_fasta_file_against_blast_database (Str $fastafilename, Str $full_path_blast_db, Alignment_Parameters $alignment_parameters) {
	my $dir = $ENV{'blast_db'};
	my $temp_file = tempdir(DIR => $dir);
	$temp_file = $temp_file.'/blastresult';
	my $factory = $alignment_parameters->create_blast_factory(8,$temp_file,$full_path_blast_db); 
	
	$factory->blastall($fastafilename);# how to specify path here?
	return $temp_file; # location of result
    } # private_blast_fasta_file_against_blast_database #

    method private_get_sequence_to_sequence_matchings_from_m8_file (Str $filename, Int $max_number_of_matches) {
	# retrieves all matching sequence-IDs, but ignores those query sequences which have
	# more than $max_number_of_matches matches
	
	my %hash = ();
	my $searchio = new Bio::SearchIO(-format => 'blasttable', -file => $filename);
	$searchio->attach_EventHandler(Bio::SearchIO::FastHitEventBuilder->new);
	while( my $result = $searchio->next_result ) {
	    if ($result->num_hits <= $max_number_of_matches) {
		my @targets;
		while( my $hit = $result->next_hit ) {
		    # hits will not have HSPs because of m8 format of BLAST output file
		    push(@targets,$hit->name);
		}
		my $query = $result->query_name;
		$GU->user_info(3,"query name: ". $result->query_name."\n");
		@{$hash{$query}} = @targets;
	    }
	    else {
	      $GU->user_info(1,"more than the desired number of BLAST records exists for this id\n");
	    }
	}
	my $result = Alignment_Set->new(matches => \%hash);
	$GU->user_info(3,"alignment set:\n");
	$GU->user_info(3,Dumper($result));
	return $result;
    } # private_get_sequence_to_sequence_matchings_from_m8_file #

} # Alignment_Utilities #
