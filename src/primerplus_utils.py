"""Collection of functions used to interface with 
   the Primer3Plus and SNPCheck websites"""

import time
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By


def primer3plus_input(browser_primers, sequence):
    """Function for interacting with Primer3Plus,
    feeding the input sequence and generating primer pairs"""

    # Input sequence
    browser_primers.get(
        "http://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi"
    )
    primer_text_form = browser_primers.find_element_by_id("sequenceTextarea")
    primer_text_form.send_keys(sequence)
    pick_primers_button = browser_primers.find_element_by_id(
        "primer3plus_pick_primers_button"
    )
    pick_primers_button.click()

    # Parse the generated primers
    forward_primers = []
    reverse_primers = []
    for i in range(10):
        if (i % 2) == 0:
            try:
                result_primer = browser_primers.find_element_by_id(
                    "".join("PRIMER_" + str(i) + "_SEQUENCE")
                )
                forward_primers.append(result_primer.get_attribute("value"))
            except:
                return "No more sets"
                # compare all sets with the least SNPs, open the link of each and check
                # the allele freqeuncy, if it's below 0.5 then download that set
        else:
            try:
                result_primer = browser_primers.find_element_by_id(
                    "".join("PRIMER_" + str(i) + "_SEQUENCE")
                )
                reverse_primers.append(result_primer.get_attribute("value"))
            except:
                pass
    product_result = browser_primers.find_elements_by_css_selector(
        "tr.primer3plus_primer_pair td"
    )
    product_sizes = [
        size.text[-6:-2] for size in product_result if "Product Size:" in size.text
    ]
    return list(zip(forward_primers, reverse_primers, product_sizes))


def SNPCheck_input(gene_name, chromosome_number, forward_primer, reverse_primer, processed_exons):
    """Function to generate the input string for SNPCheck"""

    return "".join(
        gene_name
        + "_Ex"
        + str(processed_exons[0])
        + "_F-"
        + gene_name
        + "_Ex"
        + str(processed_exons[0])
        + "_R "
        + forward_primer
        + " "
        + reverse_primer
        + " "
        + str(chromosome_number)
    )


def SNPCheck_processing(browser_SNPCheck, input_string, product_size):
    """Function to process one set of primers with SNPCheck"""

    browser_SNPCheck.get("https://genetools.org/SNPCheck/snpcheck.htm")
    SNP_text_form = browser_SNPCheck.find_element_by_id("primerPairText")
    SNP_text_form.clear()
    SNP_text_form.send_keys(input_string)
    SNP_amplicon_form = browser_SNPCheck.find_element_by_id("maxAmpliconSize")
    SNP_amplicon_form.clear()
    SNP_amplicon_form.send_keys(product_size)
    SNPCheck_button = browser_SNPCheck.find_element_by_id("snpcheckButton")
    SNPCheck_button.click()


def SNPCheck_result(browser_SNPCheck, gene_name, processed_exons):
    """Function to extract the SNP Result and amplicon size"""

    SNP_result = WebDriverWait(browser_SNPCheck, 100).until(
        EC.presence_of_element_located(
            (
                By.XPATH,
                "/html/body/div[3]/div/div/div[3]/table/tbody/tr/td/table[1]/tbody/tr[1]/td[8]/div/div[1]",
            )
        )
    )
    amplified_region_element = browser_SNPCheck.find_element_by_id(
        gene_name
        + "_Ex"
        + str(processed_exons[0])
        + "_F-"
        + gene_name
        + "_Ex"
        + str(processed_exons[0])
        + "_R.region"
    )

    amplicon = amplified_region_element.text[:3]
    print("Amplicon size =", amplicon)
    return SNP_result.text, amplicon


def specificity_check(browser_SNPCheck):
    """Function to check whether primers are specific
    in amplifying the sequence of interest only"""

    specificity_element = browser_SNPCheck.find_elements_by_css_selector(
        "tr.warnings td span"
    )
    unspecific_primer = [
        specificity.text for specificity in specificity_element]
    return unspecific_primer


def download(browser_SNPCheck):
    """Function to download the pdf and xls result files"""
    download_buttons = browser_SNPCheck.find_elements_by_class_name(
        "img_control")
    for button in download_buttons:
        if (
            button.get_attribute("src")
            == "https://genetools.org/SNPCheck/img/excel.gif"
        ):
            button.click()
        if button.get_attribute("src") == "https://genetools.org/SNPCheck/img/pdf.gif":
            button.click()
            browser_SNPCheck.find_element_by_id(
                "dijit_MenuItem_0_text").click()


def SNP_positions(SNP_element, chromosome_number):
    """Function to extract SNP positions"""

    SNP_positions = [
        SNP.text[len(chromosome_number) + 1:]
        for SNP in SNP_element
        if SNP.text.startswith(chromosome_number)
        and ":" in SNP.text
        and "KG" not in SNP.text
    ]
    for SNP in SNP_positions:
        if ".." in SNP:
            first, last = SNP.split("..")
            SNP_positions.append(first)
            SNP_positions.append(last)
    SNP_positions = [SNP for SNP in SNP_positions if ".." not in SNP]
    SNP_positions = list(dict.fromkeys(SNP_positions))
    return SNP_positions


def splitting_exon(ordered_exons, processed_exons):
    """Function to split the exon in two chunks (A & B)"""

    first_half = processed_exons[1][: len(processed_exons[1]) // 2]
    second_half = processed_exons[1][len(processed_exons[1]) // 2:]
    chunk_A = (
        first_half
        + second_half.replace("]", "")[:80]
        + "]"
        + second_half.replace("]", "")[80:300]
    )

    chunk_B = (
        first_half.replace("[", "")[-300:-80]
        + "["
        + first_half.replace("[", "")[-80:]
        + second_half
    )
    try:
        str(processed_exons[0] % 10) in "0123456789"
        ordered_exons[processed_exons[0] -
                      1] = (str(processed_exons[0]) + "a", chunk_A)
        ordered_exons = (
            ordered_exons[: processed_exons[0]]
            + [(str(processed_exons[0]) + "b", chunk_B)]
            + ordered_exons[processed_exons[0]:]
        )
    except:
        digits = [digit for digit in processed_exons[0]
                  if digit in "0123456789"]
        digits = int("".join(digits))
        letters = [
            letter
            for letter in processed_exons[0]
            if letter in "abcdefghijklmnopqrstuvwxyz"
        ]
        letters = "".join(letters)
        ele = tuple((str(digits) + letters + "a", chunk_A))
        ordered_exons[ordered_exons.index(processed_exons)] = ele
        ordered_exons = (
            ordered_exons[: ordered_exons.index(ele) + 1]
            + [tuple((str(digits) + letters + "b", chunk_B))]
            + ordered_exons[ordered_exons.index(ele) + 1:]
        )
    return ordered_exons


def last_resort(browser_SNPCheck):
    """Function to pick the best set with SNPs if all else fails"""

    SNP_links = browser_SNPCheck.find_elements_by_css_selector(
        'div[id$=".snp_table"] a'
    )
    if not SNP_links:
        print("SNPs not validated")
        return True
    else:
        for link in SNP_links:
            try:
                link.click()
                time.sleep(3)
            except:
                pass
    windows = browser_SNPCheck.window_handles
    all_freqs = []
    for tab in windows:
        if tab == windows[0]:
            continue
        browser_SNPCheck.switch_to.window(tab)
        WebDriverWait(browser_SNPCheck, 60).until(
            EC.presence_of_all_elements_located(
                (By.CSS_SELECTOR, "dl.usa-width-one-half dd div, span")
            )
        )
        frequencies = [
            frequency.text
            for frequency in browser_SNPCheck.find_elements_by_css_selector(
                "dl.usa-width-one-half dd div, span"
            )
            if ")" in frequency.text[-1:]
            and "=" in frequency.text
            or frequency.text[:4] == "None"
        ]
        all_freqs.extend(frequencies)
        try:
            browser_SNPCheck.find_element_by_id("expandfrequency").click()
            more_freqs = browser_SNPCheck.find_elements_by_css_selector(
                "div#remn_summ_freq div"
            )
            for k in more_freqs:
                all_freqs.append(k.text)
        except:
            pass
        browser_SNPCheck.close()
    accepted_freqs = 0
    for frequency in all_freqs:
        if frequency == "None":
            print("Allele frequency is:", frequency)
            accepted_freqs += 1
        elif float(frequency.replace("=", " ").split()[1]) < 0.5:
            print(frequency.replace("=", " ").split()[1], "is acceptable!")
            accepted_freqs += 1
        else:
            print("Allele frequency is too high:", frequency)
            break
    if len(all_freqs) == accepted_freqs:
        browser_SNPCheck.switch_to.window(windows[0])
        download()
        return True
    else:
        return False


# Function to pick the best set with SNPs if all else fails


def full_SNPCheck(five_primer_sets,
                  seq_to_change,
                  gene_name,
                  chromosome_number,
                  browser_SNPCheck,
                  processed_exons,
                  primer_sets,
                  ordered_full_seqs):
    """This function uses all functions defined above
    to do the following:
    - Generate primer pairs
    - Check the presence of SNPs
    - Check specificity

    If all primers have SNPs, it will try to generate better primers.
    All this will take place for an individual exon (or exon chunk).
    """

    NO_SNP_sets = []
    UNSPECIFIC_GOOD_sets = []
    filtered_out_set = []
    sets_with_SNPS = []
    emergency_sets = []
    unspecific_emergency_sets = []
    long_amplicon_sets = []
    for pair in range(len(five_primer_sets)):
        SNP_input = SNPCheck_input(
            gene_name,
            chromosome_number,
            five_primer_sets[pair][0],
            five_primer_sets[pair][1],
        )
        SNPCheck_processing(SNP_input, five_primer_sets[pair][2])
        result, amplicon_size = SNPCheck_result()

        # Scenario 1: no SNPs and immediate donwload
        if result == "No SNPs" and int(amplicon_size) <= 500:
            try:
                UNSPECIFIC_primer = specificity_check()[0].split()[8]
                UNSPECIFIC_GOOD_sets.append(SNP_input)
                total = len(specificity_check())
                if total == 1:
                    unspecific_emergency_sets.append(
                        tuple(
                            (
                                int(specificity_check()[0].split()[12]),
                                SNP_input,
                                five_primer_sets[pair][2],
                            )
                        )
                    )
                elif total == 2:
                    unspecific_emergency_sets.append(
                        tuple(
                            (
                                (
                                    int(specificity_check()[0].split()[12])
                                    + int(specificity_check()[1].split()[12])
                                ),
                                SNP_input,
                                five_primer_sets[pair][2],
                            )
                        )
                    )
                print("Good but NON-SPECIFIC Primer:", UNSPECIFIC_primer)
            except:
                print("No SNPs, downloading primer set:", SNP_input)
                NO_SNP_sets.append(SNP_input)
                download()
                break

        # Scenario 2: SNP found but irrelevant
        elif result == "All Filtered Out" and int(amplicon_size) <= 500:
            try:
                UNSPECIFIC_primer = specificity_check()[0].split()[8]
                UNSPECIFIC_GOOD_sets.append(SNP_input)
                total = len(specificity_check())
                if total == 1:
                    unspecific_emergency_sets.append(
                        tuple(
                            (
                                int(specificity_check()[0].split()[12]),
                                SNP_input,
                                five_primer_sets[pair][2],
                            )
                        )
                    )
                elif total == 2:
                    unspecific_emergency_sets.append(
                        tuple(
                            (
                                (
                                    int(specificity_check()[0].split()[12])
                                    + int(specificity_check()[1].split()[12])
                                ),
                                SNP_input,
                                five_primer_sets[pair][2],
                            )
                        )
                    )
                print("Good but NON-SPECIFIC Primer:", UNSPECIFIC_primer)
            except:
                print(
                    "Filtered out SNPs found for this primer set. Checking the next one."
                )
                filtered_out_set.append(
                    tuple((SNP_input, five_primer_sets[pair][2])))

        # Scenario 3: significant SNPs found
        elif "SNPs found in" in result and int(amplicon_size) <= 500:
            total = specificity_check()
            if len(total) == 0:
                print("SNPs FOUND!!")
                emergency_sets.append(
                    tuple((result, SNP_input, five_primer_sets[pair][2]))
                )
            else:
                print("SNPs FOUND!! Also one primer is unspecific")
            sets_with_SNPS.append(
                tuple((SNP_input, five_primer_sets[pair][2])))
        else:
            print("Error? Something wrong with this primer set:", SNP_input)
            print("Amplicon size:", amplicon_size)
            long_amplicon_sets.append(
                tuple((SNP_input, five_primer_sets[pair][2])))

    # If primers have only irrelevant SNPs, download the first set
    if filtered_out_set and not NO_SNP_sets:
        SNPCheck_processing(filtered_out_set[0][0], filtered_out_set[0][1])
        WebDriverWait(browser_SNPCheck, 100).until(
            EC.presence_of_element_located(
                (
                    By.XPATH,
                    "/html/body/div[3]/div/div/div[3]/table/tbody/tr/td/table[1]/tbody/tr[1]/td[8]/div/div[1]",
                )
            )
        )
        download()

    # If primers have no SNPs but are unspecific, try to generate an optimised result
    elif UNSPECIFIC_GOOD_sets and not NO_SNP_sets and not filtered_out_set:
        if UNSPECIFIC_primer == "1":
            Fw_sequence = list(UNSPECIFIC_GOOD_sets[0].split()[1])
            central_position = int(len(Fw_sequence) / 2)
            Fw_sequence[central_position - 1] = "N"
            seq_to_change = seq_to_change.replace(
                UNSPECIFIC_GOOD_sets[0].split()[1], "".join(Fw_sequence)
            )
            return seq_to_change, emergency_sets, unspecific_emergency_sets
        elif UNSPECIFIC_primer == "2":
            complement = {
                "A": "T",
                "a": "t",
                "T": "A",
                "t": "a",
                "C": "G",
                "c": "g",
                "G": "C",
                "g": "c",
            }
            Rv_sequence = []
            for base in UNSPECIFIC_GOOD_sets[0].split()[2][::-1]:
                Rv_sequence.append(complement[base])
            central_position = int(len(Rv_sequence) / 2)
            original_Rv_seq = "".join(Rv_sequence)
            Rv_sequence[central_position - 1] = "N"
            Rv_sequence = "".join(Rv_sequence)
            seq_to_change = seq_to_change.replace(original_Rv_seq, Rv_sequence)
            return seq_to_change, emergency_sets, unspecific_emergency_sets
        else:
            print("Unexpected condition with unspecific primers ",
                  UNSPECIFIC_GOOD_sets)

    # If there are serious SNPs, try to generate new primer sets elsewhere in the sequence
    elif sets_with_SNPS and not filtered_out_set and not NO_SNP_sets:
        SNPCheck_processing(sets_with_SNPS[0][0], sets_with_SNPS[0][1])
        WebDriverWait(browser_SNPCheck, 100).until(
            EC.presence_of_element_located(
                (
                    By.XPATH,
                    "/html/body/div[3]/div/div/div[3]/table/tbody/tr/td/table[1]/tbody/tr[1]/td[8]/div/div[1]",
                )
            )
        )
        position_element = browser_SNPCheck.find_element_by_id(
            gene_name
            + "_Ex"
            + str(processed_exons[0])
            + "_F-"
            + gene_name
            + "_Ex"
            + str(processed_exons[0])
            + "_R.region"
        )
        RV_SNP_position_element = browser_SNPCheck.find_elements_by_css_selector(
            'div[id*="_R(2).snp_table"] tr td'
        )
        FW_SNP_position_element = browser_SNPCheck.find_elements_by_css_selector(
            'div[id*="_R(1).snp_table"] tr td'
        )
        start_position = int(position_element.text.split()[2][:-2])
        end_position = int(position_element.text.split()[3])
        FW_SNP_positions = SNP_positions(FW_SNP_position_element)
        RV_SNP_positions = SNP_positions(RV_SNP_position_element)
        print("SNP position in FW primer:", FW_SNP_positions)
        print("SNP position in RV primer:", RV_SNP_positions)
        print("This is where the primer starts:", start_position)
        print("This is where the primer ends:", end_position)
        Fw_sequence = list(primer_sets[0][0])
        raw_sequence = primer_sets[0][1]
        complement = {
            "A": "T",
            "a": "t",
            "T": "A",
            "t": "a",
            "C": "G",
            "c": "g",
            "G": "C",
            "g": "c",
        }
        Rv_sequence = []
        for base in raw_sequence[::-1]:
            Rv_sequence.append(complement[base])
        if FW_SNP_positions:
            for SNPs in FW_SNP_positions:
                if (
                    end_position - int(SNPs) < len(Fw_sequence)
                    and end_position - int(SNPs) >= 0
                    or int(SNPs) - start_position < len(Fw_sequence)
                    and int(SNPs) - start_position >= 0
                ):
                    try:
                        Fw_sequence[end_position - int(SNPs)] = "N"
                    except:
                        Fw_sequence[int(SNPs) - start_position] = "N"
            Fw_sequence = "".join(Fw_sequence)
            seq_to_change = seq_to_change.replace(
                primer_sets[0][0], Fw_sequence)
        if RV_SNP_positions:
            original_Rv_seq = "".join(Rv_sequence)
            for SNPs in RV_SNP_positions:
                if (
                    int(SNPs) - start_position + 1 <= len(Rv_sequence)
                    and int(SNPs) - start_position + 1 >= 0
                    or end_position - int(SNPs) + 1 <= len(Rv_sequence)
                    and end_position - int(SNPs) + 1 >= 0
                ):
                    try:
                        Rv_sequence[-(int(SNPs) - start_position + 1)] = "N"
                    except:
                        Rv_sequence[-(end_position - int(SNPs) + 1)] = "N"
            Rv_sequence = "".join(Rv_sequence)
            seq_to_change = seq_to_change.replace(original_Rv_seq, Rv_sequence)
        return seq_to_change, emergency_sets, unspecific_emergency_sets

    # Lastly, if the amplicon was too big, try splitting the exon further
    elif len(long_amplicon_sets) == 5:
        new_ordered_full_seqs = splitting_exon(ordered_full_seqs)
        return new_ordered_full_seqs, emergency_sets, unspecific_emergency_sets

    return "Primer set obtained, next exon"
