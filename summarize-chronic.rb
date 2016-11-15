#!/usr/bin/env ruby
#
# author: scott w olesen <swo@mit.edu>

class String
  def fields
    self.chomp.split("\t")
  end
end

$input_fn = '../data/cc2011.tsv'

File.open($input_fn) do |f|
  headers = f.gets.fields
  bene = headers.shift

  # these chronic conditions get a value of 1
  single_idx = %w<AMI ALZH COPD CHF DIABETES ISCHMCHT STRKETIA>.map { |x| headers.index(x) }
  # these conditions get double points
  double_idx = %w<CHRNKIDN CNCRBRST CNCRCLRC CNCRPRST CNCRLUNG CNCRENDM>.map { |x| headers.index(x) }

  raise if single_idx.map(&:nil?).any?
  raise if double_idx.map(&:nil?).any?

  puts %w<bene comorbidity>.join("\t")

  f.each do |line|
    fields = line.chomp.split("\t")
    bene = fields.shift
    fields.map! { |x| if x.to_i == 3 then 1 else 0 end }
    val = fields.values_at(*single_idx).reduce(&:+) + 2 * fields.values_at(*double_idx).reduce(&:+)
    puts [bene, val].join("\t")
  end
end
