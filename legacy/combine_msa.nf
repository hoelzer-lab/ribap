/*Comment section: */

process combine_msa {
  label 'bioruby'
  publishDir "${params.output}/msa", mode: 'copy', pattern: "coreGenome_mafft.aln" 

  input: 
    file(aln)
    file(strain_ids)

  output:
    file("coreGenome_mafft.aln")

  script:
    """
    #!/usr/bin/env ruby

    require 'bio'

    strains = {}
    strains_file = File.open("${strain_ids}", 'r')
    strains_file.each do |line|
      if line.include?('RENAMED')
        strain = line.split(',')[1].sub('_RENAMED','').chomp
	strains[strain] = ''
      end
    end
    strains_file.close

    num = strains.keys.size

    Dir.glob("*.aln").each do |aln|
      strains_in_file = `grep ">" #{aln} | wc -l`.to_i
      next if num != strains_in_file
      Bio::FastaFormat.open(aln).each do |entry|
        id = entry.definition
	seq = entry.seq.chomp
	strains.each do |name, aln|
	  if id.start_with?(name)
	    strains[name] += seq
	    break  
	  end
	end
      end
    end

    output = File.open('coreGenome_mafft.aln','w')
    strains.each do |id, seq|
      output << ">#{id}\\n#{seq}\n"
    end
    output.close
    """
}

