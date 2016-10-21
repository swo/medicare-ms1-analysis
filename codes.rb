#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

class String
  def is_i?
    self.match(/^\d+$/)
  end
end

def ndc_triple(ndc)
  raise unless ndc.length == 11
  fields = [ndc[0..4], ndc[5..8], ndc[9..10]]
  fields.map { |x| if x.is_i? then x.to_i else x end }
end

def ndc10_to_11(package_code)
  fields = package_code.split("-")
  raise "invalid package #{package_code}" unless fields.length == 3
  fields[0].rjust(5, '0') + fields[1].rjust(4, '0') + fields[2].rjust(2, '0')
end


# make abx hash
abx = Hash.new
File.open('abx.tsv') do |f|
  headers = f.gets.chomp.split("\t")
  package_idx = headers.index('NDCPACKAGECODE')
  substance_idx = headers.index('NONPROPRIETARYNAME')

  f.each do |line|
    fields = line.chomp.split("\t")
    ndc10 = fields[package_idx]
    name = fields[substance_idx]

    # metronidazole cream exeption (maybe more)
    if ndc10 != "NA"
      ndc11 = ndc10_to_11(ndc10)
      abx[ndc11] = name
    end
  end
end


Dir.glob('data/pdesaf*.tsv').each do |fn|
  File.open(fn) do |f|
    headers = f.gets.chomp.split("\t")
    ndc_idx = headers.index('PRDSRVID')
    prescriber_idx = headers.index('CCW_PRSCRBR_ID') || headers.index('PRSCRBID')

    # swo> I don't know if these are the same
    pharm_idx = headers.index('PHARM_ID') || headers.index('PRSCRBID')

    raise "could not find prescriber in header of file #{fn}" if prescriber_idx.nil?

    f.each do |line|
      fields = line.chomp.split("\t")
      ndc = fields[ndc_idx]

      if abx.key? ndc
        puts [fn, fields[prescriber_idx], ndc, abx[ndc]].join("\t")
      end
    end
  end
end
