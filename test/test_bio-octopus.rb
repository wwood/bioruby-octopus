require 'helper'
require 'tempfile'
require 'bio-tm_hmm'

class TestBioOctopus < Test::Unit::TestCase
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

  def test_wrapper_read
    Tempfile.open('spock') do |tempfile|
      tempfile.puts ''
      tempfile.flush

      pees = Bio::Spoctopus::WrapperParser.new(tempfile.path).transmembrane_proteins
      assert_equal [], pees
    end

    Tempfile.open('spock') do |tempfile|
      tempfile.puts 'pfa|PFD0635c	I	1833	1853	outside_in'
      tempfile.flush

      pees = Bio::Spoctopus::WrapperParser.new(tempfile.path).transmembrane_proteins
      assert_equal 1, pees.length
      r = pees[0]
      assert_equal 'pfa|PFD0635c', r.name
      assert_equal 1, r.transmembrane_domains.length
      t = r.transmembrane_domains[0]
      assert_equal 1833, t.start
      assert_equal 1853, t.stop
      assert r.transmembrane_type_1?
    end

    Tempfile.open('spock') do |tempfile|
      tempfile.puts 'pfa|PFD0635c	I	1833	1853	outside_in
pfa|PFD0595c	II	2	22	inside_out
pfa|PFB0610c	No Transmembrane Domain Found
pfa|PFF1525c	Unknown	2	22	outside_in
pfa|PFF1525c	Unknown	160	180	inside_out
pfa|PFF1525c	Unknown	188	208	outside_in'
      tempfile.flush

      pees = Bio::Spoctopus::WrapperParser.new(tempfile.path).transmembrane_proteins
      assert_equal 4, pees.length
      r = pees[0]
      assert_equal 'pfa|PFD0635c', r.name
      assert_equal 1, r.transmembrane_domains.length
      t = r.transmembrane_domains[0]
      assert_equal 1833, t.start
      assert_equal 1853, t.stop
      assert r.transmembrane_type_1?

      r = pees[1]

      assert r.transmembrane_type_2?

      assert_equal 'pfa|PFB0610c', pees[2].name
      assert_equal false, pees[2].has_domain?

      assert_equal 3, pees[3].transmembrane_domains.length
    end
  end
end
