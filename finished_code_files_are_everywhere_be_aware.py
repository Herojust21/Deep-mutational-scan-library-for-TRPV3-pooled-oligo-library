#Justin Lerma 07/06/2023
# This code is to take in the finklestein generated barcodes 12 nt in length and generate a new list where the first 6 barcodes do not match for 4960 barcondes, just in case 5,000.



def barcodes_into_list():

    # tries to open a file with the barcodes and inputs them into a list in python
    try:
        with open("C:/Users/19568/PycharmProjects/pythonProject/barcodes12-1.txt", 'r') as f:
            list = [line.strip() for line in f]
        return list
    # if an error is found return this message
    except error:
        print("no file found problem in barcodes_into_lists().")
        return []





def first_six_5000(barcodes):

    prefix_set = set()  # Set to store unique prefixes of 6
    result_list = []  # List to store the final barcodes that have the first 6 nucleotides unique

    for string in barcodes:
        prefix = string[:6]
        if prefix not in prefix_set:
            prefix_set.add(prefix)
            result_list.append(string)

    return result_list





def spine_into_one_list(fastaspine):
    genes = []
    current_gene = None

    with open(fastaspine, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith('>'):
                if current_gene is not None:
                    genes.append(current_gene)

                current_gene = {'header': line[1:], 'sequence': ''}
            else:
                if current_gene is not None:
                    current_gene['sequence'] += line

        if current_gene is not None:
            genes.append(current_gene)

    return genes
    # puts all the genes from the spine output into a list for use.





def split_spineoutput_into_2_lists(genes_list, split_sequence):
    first_list = []
    second_list = []
    lost_list =[]

    for gene in genes_list:
        sequence = gene['sequence']

        # Find the index of the split sequence
        index = sequence.find(split_sequence)

        if index != -1:
            # Split the gene and create two new genes
            first_gene = {'header': gene['header'] + '_part1', 'sequence': sequence[:index + len(split_sequence)]}
            second_gene = {'header': gene['header'] + '_part2', 'sequence': sequence[index + len(split_sequence):]}

            first_list.append(first_gene)
            second_list.append(second_gene)
        else:
            # If the split sequence is not found, add the gene to the first list
            lost_list.append(gene)
    print('number of failures to split where no bsmb1 site exists: ', lost_list.append(gene))
    return first_list, second_list





def check(word, list):
    if word in list:
        print("The word is in the list!")
    else:
        print("The word is not in the list!")





def tiles_prefix_generator(list1):
    prefix_set = set()  # Set to store unique prefixes of 6
    result_list = []  # List to store the final barcodes that have the first 6 nucleotides unique

    for string in list1:
        prefix = string[:13]
        if prefix not in prefix_set:
            prefix_set.add(prefix)
            result_list.append(string)

    return result_list





def split_spineoutput_into_2_lists_but_before_cut_site(genes_list, split_sequence):
    first_list = []
    second_list = []
    lost_list = []

    for gene in genes_list:
        sequence = gene['sequence']

        # Find the index of the split sequence
        index = sequence.find(split_sequence)

        if index != -1:
            # Split the gene and create two new genes
            first_gene = {'header': gene['header'] + '_part1', 'sequence': sequence[:index]}
            second_gene = {'header': gene['header'] + '_part2', 'sequence': sequence[index + len(split_sequence):]}

            first_list.append(first_gene)
            second_list.append(second_gene)
        else:
            # If the split sequence is not found, add the gene to the first list
            lost_list.append(gene)
    print('number of failures to split where no bsmb1 site exists: ', lost_list.append(gene))
    return first_list, second_list





def check(word, list):
    if word in list:
        print("The word is in the list!")
    else:
        print("The word is not in the list!")





def maketiles1to4(list1, list2, p1, p2, bc, prefixes):
    n = 0
    barcodecounter = 0
    templist = list1

    for x in range(len(templist)+1):
        if n == 1349 or n == 2717:
            barcodecounter = 0
            if n == 2717:
                return templist
        templist[n]['sequence'] = templist[n]['sequence'] + 'ATTCC' + bc[barcodecounter] + p2 + list2[n]['sequence']
        barcodecounter += 1
        #print(templist[n]['sequence']) #checks to see if it printed it right
        n +=1





def maketiles5to6(list1, list2, p1, p2, bc, prefixes):
    n = 2717
    barcodecounter = 0
    templist = list1

    for x in range(len(templist) + 1):
        if n == 3401:
            for y in range(len(templist) + 1):
                templist[n]['sequence'] = templist[n]['sequence'] + 'ATTCC' + bc[barcodecounter] + p2 + 'c' + list2[n]['sequence']
                barcodecounter += 1
                #print(templist[n]['sequence'])  # checks to see if it printed it right
                n += 1
                if n == 4066:
                    return templist
        templist[n]['sequence'] = templist[n]['sequence'] + 'ATTCC' + bc[barcodecounter] + p2 + list2[n]['sequence']
        barcodecounter += 1
        #print(templist[n]['sequence'])  # checks to see if it printed it right
        n += 1





def maketile7(list1, list2, list3, p1, p2, bc, prefixes):
    n = 4066
    barcodecounter = 0
    templist = list1
    for x in range(len(templist) + 1):
        if n >= 4731:
            barcodecounter = 0
            return templist
        templist[n]['sequence'] = templist[n]['sequence'] + 'ATTCC' + bc[barcodecounter] + p2 + list2[n]['sequence'] + 'agatGGAGACG' + list3[n]['sequence']
        barcodecounter += 1
        #print(templist[n]['sequence'])
        #print(templist[n]['header'])# checks to see if it printed it right
        n += 1





def write_fasta_file(data_list, filename):
    with open(filename, 'w') as fasta_file:
        for data_dict in data_list:
            header = data_dict['header']
            sequence = data_dict['sequence']
            fasta_file.write(f'>{header}\n{sequence}\n')

# List of dictionaries examples
#data_list = [{'header': 'Trpv3_Mut2-36_Gln2Cys_part1', 'sequence': 'TTGGTCATGTGCTTTTCGTTGGGTGGGTACGTCTCCATTCCAACAACAACACCACTGaGAAGAGCaaGCTCTTCgaggTGCaagaagaagcgactgaagaagcgcatcttcgcggctgtgtccgagggctgcgtggaggagctgcgggaactcctacaggatctgcaggacctctgcaggaggcgccGGAGACGGCTATATCCGGGGAATCGATTCTGAGTAT'}, {'header': 'Trpv3_Mut2-36_Gln2Asp_part1', 'sequence': 'TTGGTCATGTGCTTTTCGTTGGGTGGGTACGTCTCCATTCCAACAAGAACCAGACTGaGAAGAGCaaGCTCTTCgaggGATaagaagaagcgactgaagaagcgcatcttcgcggctgtgtccgagggctgcgtggaggagctgcgggaactcctacaggatctgcaggacctctgcaggaggcgccGGAGACGGCTATATCCGGGGAATCGATTCTGAGTAT'}, {'header': 'Trpv3_Mut2-36_Gln2Ser_part1', 'sequence': 'TTGGTCATGTGCTTTTCGTTGGGTGGGTACGTCTCCATTCCAACAATACCACGACTGaGAAGAGCaaGCTCTTCgaggTCTaagaagaagcgactgaagaagcgcatcttcgcggctgtgtccgagggctgcgtggaggagctgcgggaactcctacaggatctgcaggacctctgcaggaggcgccGGAGACGGCTATATCCGGGGAATCGATTCTGAGTAT'}]

# Filename example
#filename = 'output.fasta'

# Call the function to write the FASTA file example
#write_fasta_file(data_list, filename)







def main():

    #barcode_file = "C:/Users/19568/PycharmProjects/pythonProject/barcodes12-1.txt" this is the file path that I used, notice you just need the file in the correct folder
    barcodes = barcodes_into_list() #this function takes the list of available barcodes from the finklestein paper and puts them into a list as a dictionary
    #templist = ['AACAACAACACC','AACAACAACGGT','AACAACACAAGC','AAGCCTATTCGG','AAGTAAGTCAGG'] this is an example to test the code
    print("from {} free barcodes,".format(len(barcodes)) + " you get: ")




    libraryofbc = first_six_5000(barcodes) #This function gets the first 6 unique nt from the list of the 12 nt length barcodes from the free barcodes.
    print(len(libraryofbc), "unique barcodes") # tells me how many unique barcodes there are.






    fastaspine = "C:/Users/19568/Desktop/Genesnapfiles/TRPV3/DMS library/trpv3_spine_libraryHopeVersion/All_Oligosc.txt"
    spine = spine_into_one_list(fastaspine)
    #this function splites the spineoutput into a list of dictionaries

    '''
    for s in spine:
        print('Sequence:', s['sequence'])

    print(len(spine))
    '''






    num = len([1 for line in open("C:/Users/19568/Desktop/Genesnapfiles/TRPV3/DMS library/trpv3_spine_libraryHopeVersion/All_Oligosc.txt") if line.startswith(">")])
    print('there are ', num, ' genes in this file because it does not include the wild type for each position which is 4980 minus 249.')
    print('The tiles will be split into pairs of 2s so that tiles 1-2 are a pairs 3-4 5-6 and then at last 7. \nthe wild type will be given the last barcode on the list of unique barcodes.'
          'each pair will need 1440 barcodes and 3023 divided by 2 means there are 1511.5 unique \nbarcodes to use for each tile using this method but only the first '
          '1440 will be used.')

    #the commented out code above prints the sequences of the list and the length, there should be 4980 genes but i only get 4731 because it does not give me the wild types.

    #this new section splites the spine list into two list cut by the first bsmb1 stie






    bsmb1 = 'CGTCTCC'
    bsmb1reverse = 'GGAGACG'
    padlock1 = 'ATTCC'
    padlock2 = 'ACTGaGAAGAGCaaGCTCTTC'
    counter = 0
    tiles = []
    tiles1t4 = []
    tiles5t6 = []
    tile7 = []
    genedic_handle1 =[]
    mutation_and_genetic_handle2 = []
    tile7_handle2 = []
    genedic_handle1, mutation_and_genetic_handle2 = split_spineoutput_into_2_lists(spine,bsmb1 )
    tilenames = []
    for s in genedic_handle1:
        tilenames.append(s['header'])
    #print(genedic_handle1[:3]) #wanted to inspect the list
    prefixes = tiles_prefix_generator(tilenames)
    #print(prefixes)
    tiles1t4 = maketiles1to4(genedic_handle1, mutation_and_genetic_handle2, padlock1, padlock2, libraryofbc, prefixes)
    tiles5t6 = maketiles5to6(genedic_handle1, mutation_and_genetic_handle2, padlock1, padlock2, libraryofbc, prefixes)


    #for tile 7 i must first split it before the second cut site so i can add agat for 100 percent fidelity, lucky i made a function for that
    tile7_handle21, tile7_handle22 = split_spineoutput_into_2_lists_but_before_cut_site(mutation_and_genetic_handle2, bsmb1reverse)
    #print(tile7_handle21)
    tile7 = maketile7(genedic_handle1, tile7_handle21, tile7_handle22, padlock1, padlock2, libraryofbc, prefixes)
    #tiles = create_strings(genedic_handle1, mutation_and_genetic_handle2, padlock1, padlock2, barcodes)
    #print(len(tiles))
    #print(tile7[:3]) used to see what chatgpt suggested


    #putting into fasta list
    #file1 = 'C:/Users/19568/Desktop/Genesnapfiles/TRPV3/DMS library/tiles1-4.fasta'
    file2 = 'C:/Users/19568/Desktop/Genesnapfiles/TRPV3/DMS library/tiles5-6.fasta'
    file3 = 'C:/Users/19568/Desktop/Genesnapfiles/TRPV3/DMS library/tile7.fasta'

    #fasta1to4 = write_fasta_file(tiles1t4, file1)
    fasta5to6 = write_fasta_file(tiles5t6, file2)
    fsta7 = write_fasta_file(tile7, file3)









    #print(len(genedic_handle1), len(mutation_and_genetic_handle2)) #to check if the i kept all the mutations

    '''for s in genedic_handle1:
        counter+=1
        print('Sequence:', s['sequence'])
        if counter == 3:
            counter =0
            break


    for s in mutation_and_genetic_handle2:
        counter += 1
        print('Sequence:', s['sequence'])
        if counter == 3:
            counter = 0
            break
''' # this peace of code above is to check to see if the splitting of the spine output split correctly and was eye checked by myself as you cannot test with a file yet.



    '''for s in genedic_handle1:
        counter += 1
        print('Header:', s['header'])
        if counter == 3:
            counter = 0
            break'''


if __name__ == '__main__':
    main()

