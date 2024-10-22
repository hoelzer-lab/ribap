/*Comment section: */

process generate_upsetr_input {
  label 'python3'
  publishDir "${params.output}/09-upsetr", mode: 'copy', pattern: "gene_subsets" 
  publishDir "${params.output}/09-upsetr", mode: 'copy', pattern: "*_subset.txt" 

  input:
    tuple val(ident), file(holy_table)
    file(strain_ids)

  output:
    file("gene_subsets")
    file("*_subset.txt")

  script:
"""
#!/usr/bin/env python

input_file = '${holy_table}'
input_strain_ids = '${strain_ids}'
output_file = 'gene_subsets'

with open(input_strain_ids) as file:
    strain_dict = {}
    for line in file:
        line = line.strip()
        prokka_id = line.split(',')[0]
        my_strain = line.split(',')[1]
        strain_dict[prokka_id] = my_strain

with open(input_file) as holytable:
    my_dict = {}
    my_list = []
    for line in holytable:
        line = line.strip()
        if line.startswith("Cluster_ID"):
            strain = line.split('\\t')
            my_list = strain[3:]
            for i in my_list:
                if len(${params.annotation_file}) > 1:
                    my_dict[i] = []
                else:
                    my_dict[i + '_RENAMED'] = []
        else:
            ids = line.split('\\t')
            for x in line.split('\\t')[3:]:
                if x == "NA":
                    continue
                else:
                    my_dict[strain_dict[x.split('_')[0]]].append(ids[0])
                    
for i in my_dict:
    print(i, my_dict[i])
                    
with open(output_file, 'w') as out:
    for key, value in my_dict.items():
        out.write(key + '\\n')
        for entry in value:
            out.write(entry + '\\n')

     
# write to multiple output files                
for key, value in my_dict.items():
    with open(f'{key}_subset.txt', 'w') as out:
        for entry in value:
            out.write(entry + '\\n')
    """
}

