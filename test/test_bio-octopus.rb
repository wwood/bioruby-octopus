require 'helper'
require 'tempfile'

class TestBioOctopus < Test::Unit::TestCase
	DATA_DIR = File.join(Dir.pwd,'test','data')
	 
  def test_no_tmd_result
    res = Bio::Spoctopus::Result.create_from_output([
        '>wrapperSeq',
        'gggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg',
        'ggggggggggggggggggggggggggggggggggggg'
      ].join("\n"))

    assert_kind_of Bio::Transmembrane::SignalPeptideTransmembraneDomainProtein, res
    assert_equal [], res.transmembrane_domains
    assert_equal false, res.signal?
  end

  def test_two_tmd_result
    res = Bio::Spoctopus::Result.create_from_output([
        '>wrapperSeq',
        'iiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooo',
        'ooooMMMMMMMMMMMMMMMMMMMMMiiiiiMMMMMMMMMMMMMMMMMMMMMo'
      ].join("\n"))

    assert_kind_of Bio::Transmembrane::SignalPeptideTransmembraneDomainProtein, res
    assert_equal 3, res.transmembrane_domains.length
    assert_equal 11, res.transmembrane_domains[0].start
    assert_equal 31, res.transmembrane_domains[0].stop
    assert_equal 112-1, res.transmembrane_domains[2].stop

    # test orientation
    assert_equal Bio::Transmembrane::OrientedTransmembraneDomain::INSIDE_OUT, res.transmembrane_domains[0].orientation
    assert_equal Bio::Transmembrane::OrientedTransmembraneDomain::OUTSIDE_IN, res.transmembrane_domains[1].orientation
    assert_equal Bio::Transmembrane::OrientedTransmembraneDomain::INSIDE_OUT, res.transmembrane_domains[2].orientation
  end

  def test_all_tmd_result
    res = Bio::Spoctopus::Result.create_from_output([
        '>wrapperSeq',
        'MMMMMMMMMMMMMMMMMMMMM'
      ].join("\n"))

    assert_equal 1, res.transmembrane_domains.length
    assert_equal 1, res.transmembrane_domains[0].start
    assert_equal 21, res.transmembrane_domains[0].stop
    assert_equal Bio::Transmembrane::OrientedTransmembraneDomain::UNKNOWN, res.transmembrane_domains[0].orientation
  end

  def test_tmd_at_end_result
    res = Bio::Spoctopus::Result.create_from_output([
        '>wrapperSeq',
        'oooMMMMMMMMMMMMMMMMMMMMM'
      ].join("\n"))

    assert_equal 1, res.transmembrane_domains.length
    assert_equal Bio::Transmembrane::OrientedTransmembraneDomain::OUTSIDE_IN, res.transmembrane_domains[0].orientation
  end

  def test_signal_peptide
    res = Bio::Spoctopus::Result.create_from_output([
        '>wrapperSeq',
        'nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnSSSSSSSSSSSSSSSoooooooooooooooooooooooooooooooooooooooooooooooooooo'
      ].join("\n"))

    assert res.signal?
    assert_equal false, res.has_domain?
    assert_equal 31, res.signal_peptide.start
    assert_equal 45, res.signal_peptide.stop
  end

  def test_reentrant
    res = Bio::Spoctopus::Result.create_from_output([
        '>wrapperSeq',
        'iiiirrrrrrriiiiiiiiiiiMMMMMM
MMMMMMMMMMMMMMMoooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiMMMMMMMMMMMM
MMMMMMMMMooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiMMMMMM
MMMMMMMMMMMMMMMooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiii
iiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooo'
      ].join("\n"))

    assert_equal false, res.signal?
    assert res.has_domain?
  end
  
  def test_wrapper_by_actually_running_the_underlying_program
  	sequence = "MKFASKKNNQKNSSKNDERYRELDNLVQEGNGSRLGGGSCLGKCAHVFKLIFKEIKDNIFIYILSIIYLSVCVMNKIFAK
RTLNKIGNYSFVTSETHNFICMIMFFIVYSLFGNKKGNSKERHRSFNLQFFAISMLDACSVILAFIGLTRTTGNIQSFVL
QLSIPINMFFCFLILRYRYHLYNYLGAVIIVVTIALVEMKLSFETQEENSIIFNLVLISALIPVCFSNMTREIVFKKYKI
DILRLNAMVSFFQLFTSCLILPVYTLPFLKQLHLPYNEIWTNIKNGFACLFLGRNTVVENCGLGMAKLCDDCDGAWKTFA
LFSFFNICDNLITSYIIDKFSTMTYTIVSCIQGPAIAIAYYFKFLAGDVVREPRLLDFVTLFGYLFGSIIYRVGNIILER
KKMRNEENEDSEGELTNVDSIITQ".gsub(/\n/,'') #this is PfCRT, MAL7P1.27
		Tempfile.open('pfcrt') do |tempfile|
			tempfile.puts sequence
			tempfile.close
			
			blastdb_path =  File.join(DATA_DIR,'dummyLegacyDb')
  		result = Bio::Spoctopus::Wrapper.new.calculate(sequence, blastdb_path)
			assert_kind_of Bio::Transmembrane::SignalPeptideTransmembraneDomainProtein, result
			assert_equal false, result.signal?
			assert_equal 58, result.transmembrane_domains[0].start
		end
  end
end
