# Import all modules required
import os
import time
import math
import textwrap
from selenium import webdriver
from preprocessing_utils import left_flanking, right_flanking, letter_permutations
from primerplus_utils import *

# Set the working directory and open the genetic sequence file
os.chdir("Primers Design Program")
file = open("Homo_sapiens_COL1A1_201_sequence.txt", "r")
# ------------------------------------------------------------------------------------

# STEP1: Parse the introns and exons from the genetic sequence
lowseq, upseq = "", ""
introns, exons = [], []
for line in file:
    if line[0] in "atcg":
        lowseq = lowseq + line.strip()
        if upseq == "":
            continue
        exons.append(upseq)
        upseq = ""
    elif line[0] in "ATCG":
        upseq = upseq + line.strip()
        if lowseq == "":
            continue
        introns.append(lowseq)
        lowseq = ""
introns.append(lowseq)
file.close

# STEP 2: Process exons
# Handle exons smaller than 500bp
full_seqs = []
ordering = []
for k in range(len(exons)):
    if len(exons[k]) <= 320 and len(introns[k]) >= 300 and len(introns[k + 1]) >= 300:
        full_seqs.append(
            introns[k][-300:-25]
            + "["
            + introns[k][-25:]
            + exons[k]
            + introns[k + 1][:25]
            + "]"
            + introns[k + 1][25:300]
        )
        ordering.append(k + 1)
    elif len(exons[k]) <= 320 and len(introns[k]) < 300 and len(introns[k + 1]) >= 300:
        left_seq = left_flanking(exons, introns, k)
        full_seqs.append(
            left_seq[-300: -len(introns[k])]
            + introns[k][:-25]
            + "["
            + introns[k][-25:]
            + exons[k]
            + introns[k + 1][:25]
            + "]"
            + introns[k + 1][25:300]
        )
        ordering.append(k + 1)
    elif len(exons[k]) <= 320 and len(introns[k]) >= 300 and len(introns[k + 1]) < 300:
        right_seq = right_flanking(exons, introns, k)
        full_seqs.append(
            introns[k][-300:-25]
            + "["
            + introns[k][-25:]
            + exons[k]
            + introns[k + 1][:25]
            + "]"
            + introns[k + 1][25:]
            + right_seq[len(introns[k + 1]): 300]
        )
        ordering.append(k + 1)
    elif len(exons[k]) <= 320 and len(introns[k]) < 300 and len(introns[k + 1]) < 300:
        left_seq = left_flanking(exons, introns, k)
        right_seq = right_flanking(exons, introns, k)
        full_seqs.append(
            left_seq[-300: -(len(introns[k]))]
            + introns[k][:-25]
            + "["
            + introns[k][-25:]
            + exons[k]
            + introns[k + 1][:25]
            + "]"
            + introns[k + 1][25:]
            + right_seq[len(introns[k + 1]): 300]
        )
        ordering.append(k + 1)

    # Handle exons larger than 500bp by splitting into chunks of 150bp
    else:
        parts = math.ceil(len(exons[k]) / math.ceil(len(exons[k]) / 150))
        split_exons = textwrap.wrap(exons[k], parts)
        for count, letter in enumerate(letter_permutations()):
            ordering.append("".join(str(k + 1) + letter))
            if count == len(split_exons) - 1:
                break
        for i in range(len(split_exons)):
            if i == 0:
                left_seq = left_flanking(exons, introns, k)
                right_seq = ""
                for R in split_exons[i + 1:]:
                    right_seq = right_seq + R
                right_seq = right_seq + right_flanking(exons, introns, k)
                full_seqs.append(
                    left_seq[-300:-25]
                    + "["
                    + left_seq[-25:]
                    + split_exons[i]
                    + right_seq[:70]
                    + "]"
                    + right_seq[70:300]
                )
            elif i == len(split_exons) - 1:
                right_seq = right_flanking(exons, introns, k)
                left_seq = ""
                for L in split_exons[:i]:
                    left_seq = left_seq + L
                left_seq = left_flanking(exons, introns, k) + left_seq
                full_seqs.append(
                    left_seq[-300:-70]
                    + "["
                    + left_seq[-70:]
                    + split_exons[i]
                    + right_seq[:25]
                    + "]"
                    + right_seq[25:300]
                )
            elif i != len(split_exons) - 1 and i != 0:
                left_seq = left_flanking(exons, introns, k)
                right_seq = ""
                for L in split_exons[:i]:
                    left_seq = left_seq + L
                for R in split_exons[i + 1:]:
                    right_seq = right_seq + R
                right_seq = right_seq + right_flanking(exons, introns, k)
                full_seqs.append(
                    left_seq[-300:-70]
                    + "["
                    + left_seq[-70:]
                    + split_exons[i]
                    + right_seq[:70]
                    + "]"
                    + right_seq[70:300]
                )
            else:
                print(
                    "Unusual exception when processing chunck ",
                    i + 1,
                    " of exon ",
                    k + 1,
                )
                print("Check manually")
        split_exons = []
# Combine ordering and full sequences into a list of tuples
ordered_full_seqs = list(zip(ordering, full_seqs))


# STEP 3a: Configure selenium instance profile
# Configure download settings and download path
profile = webdriver.FirefoxProfile()
profile.set_preference("browser.download.folderList", 2)
profile.set_preference(
    "browser.download.dir",
    "C:\\Users\\Alex\\Desktop\\HelloWorld\\Primers Design Program\\Primer sets",
)
profile.set_preference(
    "browser.helperApps.neverAsk.saveToDisk",
    (
        "application/vnd.ms-excel, application/pdf"
    ),  # Set MIME types for automatic download
)
profile.set_preference("pdfjs.disabled", True)  # Disable PDF.js viewer

# STEP 3b: Initiate Firefox broswer instances
browser_SNPCheck = webdriver.Firefox(firefox_profile=profile)
browser_primers = webdriver.Firefox()

# Navigate to the Ensembl page of the gene of interest
browser_primers.get(
    "http://asia.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000108821;r=17:50183289-50201632;t=ENST00000225964"
)

# Extract gene information from the Ensembl page
gene_tab = browser_primers.find_element_by_partial_link_text("Gene:")
gene_name = gene_tab.text[6:]
chromosome_tab = browser_primers.find_element_by_partial_link_text(
    "Chromosome")
if chromosome_tab.text[13] == ":":
    chromosome_number = chromosome_tab.text[11:13]
else:
    chromosome_number = chromosome_tab.text[11:12]

impossible_exons = []  # List to store exons for which primers cannot be obtained
splittings = 0  # Counter for the number of times the gene sequence is split
escape_loop = False  # Flag to escape a loop condition
current_exon = 0  # Index to iterate through ordered_full_seqs list

# Iterate through the exons
while current_exon <= len(ordered_full_seqs):
    plan_B = []
    plan_C = []
    processed_exons = ordered_full_seqs[current_exon]

    # Generate primer sets for the current exon
    primer_sets = primer3plus_input(processed_exons[1])

    # Check if primer sets are not available
    if primer_sets == "No more sets":
        print("Impossible to get primers for Exon:", processed_exons[0])
        impossible_exons.append(processed_exons[0])

    # Find primers clincal-grade for the exon
    else:
        new_seq = full_SNPCheck(primer_sets, processed_exons[1])

        # Check if the sequence is split and update ordered_full_seqs
        if type(new_seq[0]) == list and splittings < 2:
            ordered_full_seqs = new_seq[0]
            splittings += 1
            continue

        # Handle the case when primers are still not obtained after splitting twice
        elif type(new_seq[0]) == list and splittings == 2:
            # Handle the case when primers are still not obtained after splitting twice
            print("Impossible to get primers for Exon:", processed_exons[0])
            impossible_exons.append(processed_exons[0])
            ordered_full_seqs = new_seq[0]
            current_exon += 2
            splittings = 0

            # Check the escape loop condition
            if escape_loop == False:
                escape_loop = True
                continue
            else:
                current_exon += 2
                continue

        # Iterate until primer set is obtained for the current exon
        while new_seq != "Primer set obtained, next exon":
            primer_sets = primer3plus_input(new_seq[0])

            # Store alternative primer sets as a backup result (i.e. plan_B)
            plan_B.extend(new_seq[1])
            plan_B = list(dict.fromkeys(plan_B))
            plan_B.sort(key=lambda x: x[0])

            # Store extra alternative primer sets as an additional backup result (i.e. plan_C)
            plan_C.extend(new_seq[2])
            plan_C = list(dict.fromkeys(plan_C))
            plan_C.sort(key=lambda x: x[0])

            SNP_free_seq = new_seq[0]

            # Obtain a primer set from the first backup result (plan_B)
            if primer_sets == "No more sets" and plan_B:
                for pair in plan_B:
                    SNPCheck_processing(pair[1], pair[2])
                    WebDriverWait(browser_SNPCheck, 100).until(
                        EC.presence_of_element_located(
                            (
                                By.XPATH,
                                "/html/body/div[3]/div/div/div[3]/table/tbody/tr/td/table[1]/tbody/tr[1]/td[8]/div/div[1]",
                            )
                        )
                    )
                    solution = last_resort()
                    if solution:
                        break
                break

            # Obtain a primer set from the second backup result (plan_C) if plan_B failed
            elif primer_sets == "No more sets" and plan_C:
                SNPCheck_processing(plan_C[0][1], plan_C[0][2])
                WebDriverWait(browser_SNPCheck, 100).until(
                    EC.presence_of_element_located(
                        (
                            By.XPATH,
                            "/html/body/div[3]/div/div/div[3]/table/tbody/tr/td/table[1]/tbody/tr[1]/td[8]/div/div[1]",
                        )
                    )
                )
                download()
                break

            # Handle no available primers
            elif primer_sets == "No more sets" and not plan_B and not plan_C:
                print("Impossible to get primers for Exon:",
                      processed_exons[0])
                impossible_exons.append(processed_exons[0])
                break
            new_seq = full_SNPCheck(primer_sets, SNP_free_seq)
            if type(new_seq[0]) == list:
                try:
                    SNPCheck_processing(plan_C[0][1], plan_C[0][2])
                    WebDriverWait(browser_SNPCheck, 100).until(
                        EC.presence_of_element_located(
                            (
                                By.XPATH,
                                "/html/body/div[3]/div/div/div[3]/table/tbody/tr/td/table[1]/tbody/tr[1]/td[8]/div/div[1]",
                            )
                        )
                    )
                    download()
                    break
                except:
                    ordered_full_seqs = new_seq[0]
                    plan_B = []
                    plan_C = []
                    processed_exons = ordered_full_seqs[current_exon]

    # Update the current_exon index
    try:
        current_exon += 1
        splittings = 0
        escape_loop = False
        processed_exons = ordered_full_seqs[current_exon]
    except:
        break

# Print impossible exons if any
if impossible_exons:
    print("Primer3Plus cannot generate primers for Exons:")
    for i in impossible_exons:
        print(i)

# Close browser instances
browser_primers.quit()
browser_SNPCheck.quit()
