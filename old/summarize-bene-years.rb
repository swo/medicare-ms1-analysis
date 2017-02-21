#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

require 'tsv'

$benes = Hash.new { |h, k| h[k] = [] }

2011.upto(2014).each do |year|
  fn = "../data/bene_match_#{year}.txt"
  TSV[fn].each do |r|
    $benes[r['BENE_ID']] << year
  end
end

$benes.each_pair do |bene, years|
  puts ([bene] + years.sort).join("\t")
end
