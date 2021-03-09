import os
import pandas as pd
import random
from Bio.Seq import MutableSeq
from Bio.Data.CodonTable import unambiguous_dna_by_id
from scipy.stats import spearmanr
from scipy.stats import pearsonr
from collections import Counter
from datetime import datetime

start_time = datetime.now()

number_of_generations = 150     # Количество поколений

species_list = ['Mus_musculus']     # Список видов. Скрипт проработает про всем
                                    # 'Abbottina_rivularis', 'Uroplatus_ebenaui', 'Eopsaltria_australis',
                                    # 'Tylototriton_verrucosus', 'Cavia_porcellus' - то, что проверял


codon_table = 2     # Вариации генетического кода с NCBI:
                    # 1: Standard   2: Vertebrate Mitochondrial 3: Yeast Mitochondrial
                    # 4: Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma
                    # 5: Invertebrate Mitochondrial     6: Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
                    # 9: Echinoderm Mitochondrial; Flatworm Mitochondrial   10: Euplotid Nuclear
                    # 11: Bacterial, Archaeal and Plant Plastid     12: Alternative Yeast Nuclear
                    # 13: Ascidian Mitochondrial    14: Alternative Flatworm Mitochondrial
                    # 16: Chlorophycean Mitochondrial   21: Trematode Mitochondrial
                    # 22: Scenedesmus obliquus Mitochondrial Code   23: Thraustochytrium mitochondrial code
                    # 24: Pterobranchia Mitochondrial   25: Candidate Division SR1; Gracilibacteria
                    # 26: Pachysolen tannophilus Nuclear Code   27: Karyorelict Nuclear     28: Condylostoma Nuclear
                    # 29: Mesodinium Nuclear    30: Peritrich Nuclear  31: Blastocrithidia Nuclear

codons64 = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC',
             'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT',
             'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC',
             'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
             'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']    # Список всех возможных кодонов

for p in range(len(species_list)): # Цикл по списку видов
    start_time1 = datetime.now()
    species = species_list[p] # Берет один вид из списка
    MutSpecData = pd.read_csv('VertebratePolymorphisms.MutSpecData.OnlyFourFoldDegAllGenes.csv', sep=' ')
    probability = 0
    mutspec = list(MutSpecData[MutSpecData.Species == species].iloc[0].values[1:])  # Вытаскивает мутспек для текущего вида
    mutspec_stack = []

    for i in range(len(mutspec)):   # Представляет мутспек в виде отрезка от 0 до 1, чтобы работала рулетка
        mutspec_stack.append(probability)
        probability += float(mutspec[i])


    for sample in range(1): # Количество генераций каждого вида(см. сколько папок уже есть в папке каждого вида)
        stage_time = datetime.now()
        if sample % 2 == 0:  # Чекпоинты прогресса скрипта
            print('Цикл {} для вида {} начался  {}'.format(sample, species, datetime.now( )))
        os.makedirs("Results_Bio/{}/{}".format(species, sample), exist_ok=True)  # Создает папку для результатов

        # Генератор генома на основе CodonUsage текущего вида

        AA = pd.read_csv('AminoAcidsUsage.csv', sep='\t')
        DNA = []

        aminoacid_codon = {}     # Создает словарь {Аминокислота : [кодоны]} на основе таблиц кодонов с NCBI
        for c in unambiguous_dna_by_id[2].forward_table:
            if unambiguous_dna_by_id[2].forward_table[c] not in aminoacid_codon:
                aminoacid_codon[unambiguous_dna_by_id[2].forward_table[c]] = []
            aminoacid_codon[unambiguous_dna_by_id[2].forward_table[c]].append(c)
        aminoacid_codon['_'] = unambiguous_dna_by_id[2].stop_codons

        for l in range(10):     # Десятикратно увеличивает геном для уменьшения "стохастики" (думаю переделать:
                                # не увеличивать геном в 10 раз, а увеличить количество прогонов скрипта)
            for aminoacid in AA.columns.values[:-1]:
                for k in range(AA[AA.Species == species][aminoacid].values[0]):  # Цикл длиной кол-ву кодонов для аминокислоты
                    DNA.append(random.choice(aminoacid_codon[aminoacid]))   # Выбирает случайный кодон аминокислоты и добавляем в DNA

        CodonUsage = pd.DataFrame(columns=codons64)     # Создаем датафреймы для подсчета CodonUsage и
        Frequency = pd.DataFrame(columns=codons64)
        observed_frequency = pd.read_csv('CodonFrequency.csv', sep='\t')
        Frequency.loc[0] = list(observed_frequency[observed_frequency.Species == species].values[0][1:])
            # Записывает частоты кодонов из живого генома текущего вида

        mutations = ['A_T', 'A_G', 'A_C', 'T_A', 'T_G', 'T_C', 'G_A', 'G_T', 'G_C', 'C_A', 'C_T', 'C_G']
            # Список однонуклеотидных замен

        for generation in range(number_of_generations):
            # Мутагенез
            for codon_id in range(len(DNA)):        # Цикл по всей ДНК
                codon = MutableSeq(DNA[codon_id])   # Выбираем отдельный кодон
                mut_id = random.randint(0, 2)       # Выбираем нуклеотид из этого кодона
                mut_probability = random.uniform(0, mutspec_stack[-1])  # Задаем значение рулетки
                counter = 0
                while mut_probability >= mutspec_stack[counter + 1]:    # Смотрим на какою именно замену она выпала
                    counter += 1
                if codon[mut_id] == mutations[counter][0]:              # Сопоставляем нуклеотид замены из мутспека с
                    new_codon = codon                                   # выбранным нуклеотидом из кодона
                    codon = codon.toseq()
                    new_codon[mut_id] = mutations[counter][2]           # Подставляем нуклеотид в кодоне
                    if codon.translate(table=codon_table) == new_codon.toseq().translate(table=codon_table):
                        DNA[codon_id] = str(new_codon)      # Проверка на синонимичность мутировавшего кодона

            if generation % 1 == 0: # Регулирует раз в сколько поколений считать частоты кодонов
                # Разбивка на список кодонов
                # codons_of_seq = []
                # for l in range(len(seq)//3):
                #     codons_of_seq.append(seq.toseq()[l*3:l*3+3])
                # print('Время разбивки: {}'.format(datetime.now( ) - stage_time))
                codon_usage = Counter(DNA)      # Расчет CodonUsage
                if len(list(codon_usage)) < 64:                     # Проверка наличия всех 64 кодонов
                    for t in codons64:                              # без нее возникали ошибки
                        if t not in list(codon_usage):              # если кодон исчезает приписываем ему ноль
                            codon_usage[t] = 0
                for i in range(64):     # Записывает CodonUsage
                    CodonUsage.at[generation, CodonUsage.columns[i]] = codon_usage[CodonUsage.columns[i]]
                frequency_values = []  # Список значений частот кодонов
                total = len(DNA)   # Сумма всех кодонов
                for e in range(64):
                    frequency_values.append(CodonUsage.loc[generation].values[e] / total)   # Расчитывает частоту каждого кодона
                Frequency.loc[generation + 1] = frequency_values  # Записывает в строку текущего поколения
                CodonUsage.to_csv('Results_Bio\{}\{}\Codons.csv'.format(species, sample), sep='\t', index=False)
                Frequency.to_csv('Results_Bio\{}\{}\Frequencies.csv'.format(species, sample), sep ='\t', index=False)
                # Сохраняет значения

        # Расчет корреляций и p-value
        data = pd.read_csv('Results_Bio\{}\{}\Frequencies.csv'.format(species, sample), sep='\t')
        spearman = []
        pearson = []
        pv_s = []
        pv_p = []
        for i in data.index.values:
            data1, data2 = data.loc[0], data.loc[i]
            corrp, p_p = pearsonr(data1, data2)
            corrs, p_s = spearmanr(data1, data2)
            spearman.append(corrs)
            pearson.append(corrp)
            pv_s.append(p_s)
            pv_p.append(p_p)
        spear_corr = pd.DataFrame(index=[i for i in range(number_of_generations+1)], columns=['spearman', 'pv_s', 'pearson', 'pv_p'])
        spear_corr['spearman'] = spearman
        spear_corr['pearson'] = pearson
        spear_corr['pv_s'] = pv_s
        spear_corr['pv_p'] = pv_p
        spear_corr.to_csv('Results_Bio\{}\{}\Spearman.csv'.format(species, sample), sep='\t', index=False)    # Расчет корреляций закончен

    print('Расчет {} закончен'.format(species))
    print(datetime.now() - start_time1)


    # Расчет средних частот кодонов для ВСЕХ смодулированных образцов и корреляций
    samples_number = len(os.listdir('Results_Bio\{}'.format(species)))-2    # Количество смодулированных образцов

    species = species_list[p]
    # Средние частоты
    FrequenciesMean = pd.read_csv('Results_Bio\{}\{}\Frequencies.csv'.format(species, 0), sep='\t').copy()
        # Записывает частоты живого генома текущего вида в датафрейм текущего вида

    for sample in range(1, samples_number):
        freq = pd.read_csv('Results_Bio\{}\{}\Frequencies.csv'.format(species, sample), sep='\t').copy()
        for i in codons64:
            FrequenciesMean[i] += freq[i]
    for i in codons64:
        FrequenciesMean[i] = FrequenciesMean[i] / samples_number
    FrequenciesMean.to_csv('Results_Bio\{}\FrequenciesMean.csv'.format(species), sep='\t', index=False)

    # Корреляции
    data = pd.read_csv('Results_Bio\{}\FrequenciesMean.csv'.format(species), sep='\t').copy()
    spearman = []
    pearson = []
    pv_s = []
    pv_p = []
    for i in data.index.values[1:]:
        data1, data2 = data.loc[0], data.loc[i]
        corrp, p_p = pearsonr(data1, data2)
        corrs, p_s = spearmanr(data1, data2)
        spearman.append(corrs)
        pearson.append(corrp)
        pv_s.append(p_s)
        pv_p.append(p_p)
    spear_corr = pd.DataFrame(index=[i for i in range(number_of_generations)]
                              , columns=['spearman', 'pv_s', 'pearson', 'pv_p'])
    spear_corr['spearman'] = spearman
    spear_corr['pearson'] = pearson
    spear_corr['pv_s'] = pv_s
    spear_corr['pv_p'] = pv_p
    spear_corr.to_csv('Results_Bio\{}\SpearmanMean.csv'.format(species), sep='\t')

print('Расчет средних корреляций закончен')
print(datetime.now() - start_time)
