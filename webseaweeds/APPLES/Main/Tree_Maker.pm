### (c) copyright University of Warwick 2009, 2010 ###
### This file is part of the APPLES (Analysis of Plant Promoter-Linked Elements) package. ###

### Tree_Maker Class ###

### provide constructors for commonly used trees (use stuff from make_tree.pl)

use MooseX::Declare;

class Tree_Maker {

} # Tree_Maker #

class Evolutionary_Tree_Maker extends Tree_Maker {
    use APPLES_Datatypes qw (EvolutionaryTreeKingdom);
    
    method make_evolutionary_tree(EvolutionaryTreeKingdom $kingdom) {
	
	if ($kingdom eq 'plants') {
	    my @empty_array;

	    my $a = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'arabidopsis' );

	    my $b = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'brassica' );

	    my $c = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'papaya' );

	    my $d = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'cotton' );

	    my $e = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'medicago' );

	    my $f = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'soybean' );

	    my $g = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'castor' );

	    my $h = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'poplar' );

	    my $i = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'grape' );

	    my $j = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'tomato' );

	    my $k = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'rice' );

	    my $l = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'indica' );

	    my $m = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'glaberrima' );

	    my $n = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'sorghum' );

	    my $o = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'maize' );

	    my @ab = ($a,$b);
	    my $ab = Evolutionary_Tree->new( sub_trees => \@ab);

	    my @ab_c = ($ab, $c);
	    my $abc = Evolutionary_Tree->new( sub_trees => \@ab_c );

	    my @abc_d = ($abc, $d);
	    my $abcd = Evolutionary_Tree->new( sub_trees => \@abc_d );

	    my @ef = ($e, $f);
	    my $ef = Evolutionary_Tree->new( sub_trees => \@ef );

	    my @abcd_ef = ($abcd, $ef);
	    my $abcdef = Evolutionary_Tree->new( sub_trees => \@abcd_ef );

	    my @gh = ($g, $h);
	    my $gh = Evolutionary_Tree->new ( sub_trees => \@gh );

	    my @abcdef_gh = ($abcdef, $gh);
	    my $abcdefgh = Evolutionary_Tree->new ( sub_trees => \@abcdef_gh );

	    my @abcdefgh_i = ($abcdefgh, $i);
	    my $abcdefghi = Evolutionary_Tree->new ( sub_trees => \@abcdefgh_i );

	    my @abcdefghi_j = ($abcdefghi, $j);
	    my $abcdefghij = Evolutionary_Tree->new ( sub_trees => \@abcdefghi_j );

	    my @klm = ($k, $l, $m);
	    my $klm = Evolutionary_Tree->new ( sub_trees => \@klm );

	    my @no = ($n, $o);
	    my $no = Evolutionary_Tree->new ( sub_trees => \@no );

	    my @klm_no = ($klm, $no);
	    my $klmno = Evolutionary_Tree->new ( sub_trees => \@klm_no );

	    my @abcdefghij_klmno = ($abcdefghij, $klmno);
	    my $plant_tree = Evolutionary_Tree->new ( sub_trees => \@abcdefghij_klmno );
	}
	
	elsif  ($kingdom eq 'fungi') {
	    die 'cannot build an evolutionary tree for that kingdom yet!'; 
	}

	elsif  ($kingdom eq 'vertebrates') {
	    use APPLES_Datatypes qw(APPLESSpeciesName);
	    use constant {FALSE => 0,
			  TRUE	=> 1};
	    my @species = ('human', 'chimpanzee', 'gorilla', 'orangutan', 'macaque', 'baboon', 'colobus', 'tarier', 'night monkey', 'marmoset', 'dusky titi monkey', 'mouse lemur', 'bushbaby',
			   'tree shrew', 'rabbit', 'pika', 'squirrel', 'mouse', 'rat', 'kangaroo rat', 'guinea pig', 'cow', 'dolphin', 'alpaca', 'pig', 'dog', 'cat', 'horse', 'microbat', 'megabat',
			   'hedgehog', 'shrew', 'hyrax', 'elephant', 'lesser hedgehog', 'sloth', 'armadillo', 'wallaby', 'opossum', 'platypus', 'chicken', 'quail', 'zebrafinch', 'lizard', 'xenopus',
			   'tetraodon', 'fugu', 'medaka', 'stickleback', 'zebrafish', 'savignyi', 'intestinalis', 'lamprey');
	    # establish all species leaves for species in the list
	    my @empty_array;
	    my $a = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'human' );
	    my $b = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'chimpanzee' );
	    my $c = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'gorilla' );
	    my $d = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'orangutan' );
	    my $e = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'macaque' );
	    my $f = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'baboon' );
	    my $g = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'colobus' );
	    my $h = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'tarier' );
	    my $i = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'night monkey' );
	    my $j = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'marmoset' );
	    my $k = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'dusky titi monkey' );
	    my $l = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'mouse lemur' );
	    my $m = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'bushbaby' );
	    my $n = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'tree shrew' );
	    my $o = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'rabbit' );
	    my $p = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'pika' );
	    my $q = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'squirrel' );
	    my $r = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'mouse' );
	    my $s = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'rat' );
	    my $t = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'kangaroo rat' );
	    my $u = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'guinea pig' );
	    my $v = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'cow' );
	    my $w = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'dolphin' );
	    my $x = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'alpaca' );
	    my $y = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'pig' );
	    my $z = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'dog' );
	    my $aa = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'cat' );
	    my $bb = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'horse' );
	    my $cc = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'microbat' );
	    my $dd = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'megabat' );
	    my $ee = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'hedgehog' );
	    my $ff = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'shrew' );
	    my $gg = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'hyrax' );
	    my $hh = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'elephant' );
	    my $ii = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'lesser hedgehog' );
	    my $jj = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'sloth' );
	    my $kk = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'armadillo' );
	    my $ll = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'wallaby' );
	    my $mm = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'opossum' );
	    my $nn = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'platypus' );
	    my $oo = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'chicken' );
	    my $pp = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'quail' ); 
	    my $qq = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'zebrafinch' );
	    my $rr = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'lizard' );
	    my $ss = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'xenopus' );
	    my $tt = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'salamander' );
	    my $uu = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'tetraodon' );
	    my $vv = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'fugu' );
	    my $ww = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'medaka' );
	    my $xx = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'stickleback' );
	    my $yy = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'zebrafish' );
	    my $zz = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'lamprey' ); 
	    my $aaa = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'savignyi' );
	    my $bbb = Evolutionary_Tree->new( sub_trees => \@empty_array, root => 'intestinalis' );

	    # create subtrees
	    my @a_b = ($a, $b);
	    my $ab = Evolutionary_Tree->new( sub_trees => \@a_b );
	    my @ab_c = ($ab, $c);
	    my $abc = Evolutionary_Tree->new( sub_trees => \@ab_c );

	    my @hominids = ($abc, $d);
	    my $hominids = Evolutionary_Tree->new( sub_trees => \@hominids );

	    my @ef = ($e, $f);
	    my $ef = Evolutionary_Tree->new( sub_trees => \@ef );

	    my @old_world_monkeys = ($ef, $g);
	    my $old_world_monkeys = Evolutionary_Tree->new( sub_trees => \@old_world_monkeys );

	    my @catarrhini = ($hominids, $old_world_monkeys);
	    my $catarrhini = Evolutionary_Tree->new( sub_trees => \@catarrhini );

	    my @p_2 = ($catarrhini, $h); # tarsier
	    my $p_2 = Evolutionary_Tree->new( sub_trees => \@p_2 );

	    my @ij = ($i, $j);
	    my $ij = Evolutionary_Tree->new( sub_trees => \@ij );

	    my @new_world_monkeys = ($ij, $k);
	    my $new_world_monkeys = Evolutionary_Tree->new( sub_trees => \@new_world_monkeys );

	    my @p_1 = ($p_2, $new_world_monkeys);
	    my $p_1 = Evolutionary_Tree->new( sub_trees => \@p_1 );

	    my @strepsirrhini = ($l, $m);
	    my $strepsirrhini = Evolutionary_Tree->new( sub_trees => \@strepsirrhini );

	    my @primates = ($p_1, $strepsirrhini);
	    my $primates = Evolutionary_Tree->new( sub_trees => \@primates );

	    my @lagomorpha = ($o, $p);
	    my $lagomorpha = Evolutionary_Tree->new( sub_trees => \@lagomorpha );

	    my @muroidea = ($r, $s);
	    my $muroidea = Evolutionary_Tree->new( sub_trees => \@muroidea );

	    my @r_1 = ($q, $muroidea, $t);
	    my $r_1 = Evolutionary_Tree->new( sub_trees => \@r_1 );

	    my @rodents = ($u, $r_1);
	    my $rodents = Evolutionary_Tree->new( sub_trees => \@rodents );

	    my @m_3a = ($n, $lagomorpha, $rodents);
	    my $m_3a = Evolutionary_Tree->new( sub_trees => \@m_3a );

	    my @m_3 = ($m_3a, $primates);
	    my $m_3 = Evolutionary_Tree->new( sub_trees => \@m_3 );

	    my @m_2ad = ($v, $w, $x);
	    my $m_2ad = Evolutionary_Tree->new( sub_trees => \@m_2ad );

	    my @m_2ac = ($m_2ad, $y);
	    my $m_2ac = Evolutionary_Tree->new( sub_trees => \@m_2ac );

	    my @carnivora = ($z, $aa);
	    my $carnivora = Evolutionary_Tree->new( sub_trees => \@carnivora );

	    my @m_2ab = ($m_2ac, $carnivora);
	    my $m_2ab = Evolutionary_Tree->new( sub_trees => \@m_2ab );

	    my @chiroptera = ($cc, $dd);
	    my $chiroptera = Evolutionary_Tree->new( sub_trees => \@chiroptera );

	    my @m_2aa = ($bb, $chiroptera, $m_2ab );
	    my $m_2aa = Evolutionary_Tree->new( sub_trees => \@m_2aa );

	    my @insectivora = ($ee, $ff);
	    my $insectivora = Evolutionary_Tree->new( sub_trees => \@insectivora );

	    my @m_2a = ($insectivora, $m_2aa);
	    my $m_2a = Evolutionary_Tree->new( sub_trees => \@m_2a );

	    my @m_2b = ($gg, $hh, $ii);
	    my $m_2b = Evolutionary_Tree->new( sub_trees => \@m_2b );

	    my @edentata = ($jj, $kk);
	    my $edentata = Evolutionary_Tree->new( sub_trees => \@edentata );

	    my @m_2 = ($m_2a, $m_2b, $edentata, $m_3);
	    my $m_2 = Evolutionary_Tree->new( sub_trees => \@m_2 );

	    my @marsupials = ($ll, $mm);
	    my $marsupials = Evolutionary_Tree->new( sub_trees => \@marsupials );

	    my @m_1 = ($m_2, $marsupials);
	    my $m_1 = Evolutionary_Tree->new( sub_trees => \@m_1 );

	    my @mammals = ($m_1, $nn);
	    my $mammals = Evolutionary_Tree->new( sub_trees => \@mammals );

	    my @birds_1 = ($oo, $pp);
	    my $birds_1 = Evolutionary_Tree->new( sub_trees => \@birds_1 );

	    my @birds = ($birds_1, $qq);
	    my $birds = Evolutionary_Tree->new( sub_trees => \@birds );

	    my @diapsida = ($birds, $rr);
	    my $diapsida = Evolutionary_Tree->new( sub_trees => \@diapsida );

	    my @amniota = ($diapsida, $mammals);
	    my $amniota = Evolutionary_Tree->new( sub_trees => \@amniota );

	    my @amphibians = ($ss, $tt);
	    my $amphibians = Evolutionary_Tree->new( sub_trees => \@amphibians );

	    my @terrestrial_vertebrates = ($amniota, $amphibians);
	    my $terrestrial_vertebrates = Evolutionary_Tree->new( sub_trees => \@terrestrial_vertebrates );

	    my @fish_1a = ($uu, $vv);
	    my $fish_1a = Evolutionary_Tree->new( sub_trees => \@fish_1a );

	    my @fish_1b = ($ww, $xx);
	    my $fish_1b = Evolutionary_Tree->new( sub_trees => \@fish_1b );

	    my @fish_1 = ($fish_1a, $fish_1b);
	    my $fish_1 = Evolutionary_Tree->new( sub_trees => \@fish_1 );

	    my @fish = ($fish_1, $yy);
	    my $fish = Evolutionary_Tree->new( sub_trees => \@fish );

	    my @vertebrata = ($terrestrial_vertebrates, $fish);
	    my $vertebrata = Evolutionary_Tree->new( sub_trees => \@vertebrata );

	    my @craniata = ($vertebrata, $zz);
	    my $craniata = Evolutionary_Tree->new( sub_trees => \@craniata );

	    my @urochordata = ($aaa, $bbb);
	    my $urochordata = Evolutionary_Tree->new( sub_trees => \@urochordata );

	    my @chordata = ($craniata, $urochordata);
	    my $chordata = Evolutionary_Tree->new( sub_trees => \@chordata );

	    return $chordata;
	}

	elsif  ($kingdom eq 'insects') {
	   die 'cannot build an evolutionary tree for that kingdom yet!'; 
	}
	
	elsif  ($kingdom eq 'all') {
	   die 'cannot build an evolutionary tree for that kingdom yet!'; 
	}
	
	else {
	    die 'cannot build an evolutionary tree for that kingdom yet!';
	}
    } # make_evolutionary_tree #

} # Evolutionary_Tree_Maker #

