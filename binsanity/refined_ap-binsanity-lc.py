def refined_ap(array,names,file_name,damping,iterations,convergence,preference,path,output_directory,prefix):
    """Uses affinity propagation to make putative bins"""    
    if os.path.isdir(os.path.join(str(output_directory),str(prefix)+"-REFINED-BINS")) is False:
        os.mkdir(os.path.join(str(output_directory),str(prefix)+"-REFINED-BINS"))
    name_of_output_file = os.path.join(str(output_directory),str(prefix)+"-REFINED-BINS")
    apclust = AffinityPropagation(damping=float(damping), max_iter=int(iterations), convergence_iter=int(convergence), copy=True, preference=int(preference), affinity='euclidean', verbose=False).fit_predict(array)
    outfile_data = {}
    i = 0
    while i < len(names):
        if apclust[i] in outfile_data.keys():
            outfile_data[apclust[i]].append(names[i])
        if apclust[i] not in outfile_data.keys():
            outfile_data[apclust[i]] = [names[i]]
        i += 1
    out_name = file_name.rsplit(".",1)[0]
    with open(os.path.join(path,file_name),"r") as input2_file: 
        fasta_dict = create_fasta_dict(input2_file)                
        count = 0                
        for k in outfile_data:
            if len(outfile_data[k]) >= 5:
                output_file = open(os.path.join(name_of_output_file,str(out_name)+"-refined_%s.fna" % (k)), "w" )
                for x in outfile_data[k]:
                    output_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
                output_file.close()
                count = count + 1
            elif len(outfile_data[k]) < 5:
                if any((len(fasta_dict[x])>50000) for x in outfile_data[k]):
                    output_file = open(os.path.join(name_of_output_file,str(out_name)+"-refined_%s.fna" % (k)), "w" )
                    for x in outfile_data[k]:
                        output_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
                    output_file.close()
                    count = count +1
            else:
                output_file=open(os.path.join(output_directory,'unclustered.fna'),'a')
                for x in outfile_data[k]:
                    output_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
                output_file.close()
            if os.path.isfile(os.path.join(output_directory,'unclustered.fna')) is True:
                contig_number = 0
                for record in SeqIO.parse(os.path.join(output_directory,'unclustered.fna'),"fasta"):
                        contig_number +=1
                if contig_number < 2:
                        os.remove(os.path.join(output_directory,'unclustered.fna'))
            print "          Cluster "+str(k)+": "+str(len(outfile_data[k]))
        print ("""          Total Number of Bins: %i""" % count)
