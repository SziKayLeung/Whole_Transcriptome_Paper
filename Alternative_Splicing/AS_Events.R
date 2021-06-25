#!/usr/bin/env Rscript
# Szi Kay Leung: sl693@exeter.ac.uk
# Aim: Create functions for processing alternative splicing events in Iso-Seq: 
	# AF_AL -> Alternative First and Alternative Last Exons 
	# A5_A3 -> Alternative 5' and 3' 


#************************************* FUNCTIONS

# AF_AL <input_dat>
# Aim: Detect alternative first and alternative last exons by ordering the coordinates of isoforms (allow wobble of 10bp for collapse)
# input_dat = gtf file with coordinate of the exons and transcripts of PB.ID isoforms (SQANTI2 filtered gtf file)
AF_AL <- function(input_dat){
  # first and last exon differs between strands 
  
  ## plus 
  # exon1|--->---intron--->-----|exon2|-->---- (coordinates 1-2-3)
  # arrange coordinates in ascending order from 1, 2, 3 i.e 1st row would be 1st exon, last row would be last exon
  
  ## minus 
  # exon2|---<---intron----<----|exon1|---<--- (coordinates 1-2-3)
  # arrange coordinates in descending order from 3,2,1 i.e 1st row would be 1st exon, last row would be last exon
  
  plus <- input_dat %>% filter(source == "exon", strand == "+") %>% group_by(isoform) %>% arrange(start, .by_group = T)
  minus <- input_dat %>% filter(source == "exon", strand == "-") %>% group_by(isoform) %>% arrange(-start, .by_group = T)
  
  first_exon <- function(strand_input){strand_input %>% filter(row_number()==1)}
  last_exon <- function(strand_input){strand_input  %>% filter(row_number()==n())}
  
  ## AF
  # bind the first exon of the plus and minus strand, then group by the gtf coordinates of the first exon per gene 
  # some isoforms of the same gene would have the same coordinates 
  first_exon_output <- rbind(first_exon(plus), first_exon(minus)) 

  ## AL
  # similar to AF but only capture last exon 
  last_exon_output <- rbind(last_exon(plus), last_exon(minus)) 
  
  # Filter for wobbles in the first exon at the 5' end 
  #1. |exon1|---->-
  #2. |exon1|---->-
  #3.    |exon 1|-->----
  #4.               |exon 1|-->---
  # if AF, no difference between "start" coordinates of 1st exon after sorting (i.e. isoform 4 and isoform 1)
  # However, false positives from wobble i.e. isoform 3 and isoform 1, therfore absolute difference between start coordinates >10bp
  # Criteria based off from TAMA collapse default setting with wobble
  # Absolute difference to account for minus and positive strands after sorting
  # Remove NA as this means only one isoform for that gene or all the isoforms have the same first exon, therefore no AF
  AF <- first_exon_output %>% group_by(associated_gene) %>% arrange(start, .by_group = TRUE) %>% 
    mutate(Diff = abs(start - lag(start))) %>% filter(!is.na(Diff) & Diff > 10) %>%
    mutate(Event = "AF")
  
  # filter for wobbles in the last exon at the 3'end
  # similar concept as 5' end, but use the "end" coordinates of last exon after sorting
  #1. ----->-------|exon20|
  #2. ----->-------|exon20|
  #3. ----->  ----|exon20|
  #4. |exon20|
  # Allow wobble with 10bp between the end
  AL <- last_exon_output %>% group_by(associated_gene) %>% arrange(end, .by_group = TRUE) %>% mutate(Diff = abs(end - lag(end))) %>% 
  filter(!is.na(Diff) & Diff > 10) %>% mutate(Event = "AL")
  
  AF_AL <- rbind(AF,AL)
  return(AF_AL)
}

# A5_A3 <input_dat> <threshold>
# Aim: Detect alternative 5' and alternative 3' exons by ordering the coordinates of isoforms (allow user-defined threshold for collapse)
# input_dat = gtf file with coordinate of the exons and transcripts of PB.ID isoforms (SQANTI2 filtered gtf file)
A5_A3 <- function(input_dat, threshold){
  
  # Method 1 
  # plus strand = A5'
  #1.NO         (a*)|exon1|(a)---->------(d*)|exon 2|(d)----
  #2.NO         (b*)|exon1|(b)---->------(e*)|exon 2|(e)----
  #3.YES        (c*)|exon1   |(c)---->---(f*)|exon 2|(f)----
  
  # minus strand = A3'
  #1.NO        (a*)|exon3|(a)----<------(e*)|exon 2|(e)----
  #2.NO        (b*)|exon3|(b)----<------(f*)|exon 2|(f)----
  #3.YES       (c*)|exon3   |(c)----<---(g*)|exon 2|(g)----
  #4.NO        (d*)|exon3   |(d)----<---------------------(h*)|exon 2|(h)----
  
  
  # Method 2 
  # plus strand = A3'
  #1.NO         (a*)|exon1|(a)---->------(d*)|exon 2     |(d)----
  #2.NO         (b*)|exon1|(b)---->------(e*)|exon 2     |(e)----
  #3.YES        (c*)|exon1|(c)---->-----------(f*)|exon 2|(f)----
  
  # minus strand = A5'
  #1.NO        (a*)|exon3|(a)----<------(e*)|exon 2     |(e)----
  #2.NO        (b*)|exon3|(b)----<------(f*)|exon 2     |(f)----
  #3.YES       (c*)|exon3|(c)----<-----------(g*)|exon 2|(g)----
  #4.NO        (d*)|exon3|(d)----<------------------------(h*)|exon 2|(h)----
  
  
  # Note Isoform 4 may be pick up as A5'/A3' however it's actually just the next exon along, so will have a massive difference
  # However, need to filter out, therefore select threshold
  
  # By sorting out the end/start coordinates from ascending order: (in plus strand)
  # (a*)|exon1|(a)
  # (b*)|exon1|(b)
  # (c*)|exon1   |(c)
  #                           (d*)|exon 2|(d)
  #                           (e*)|exon 2|(e)
  #                           (f*)|exon 2|(f)
  #                           (g*)|exon 2|(g) 
  
  
  # For Method 1, A5' or A3' have to differ between the end coordinates of the same exon, but the same for start coordinates the next exon along of isoform, and same start coordinates of the same exon
  # i.e Isoform 3 and Isoform 1, whereby abs(c-b) != 0 but abs(g*-f*) == 0
  # To calculate difference between end coordinates (c-b) => lag the end coordinates
  # To calculate difference between start coordinates for the next along of isoform:
  # by isoform, create a new column for the next start coordinates
  input_dat <- input_dat %>% group_by(isoform) %>% filter(source == "exon") %>% arrange(start, .by_group = TRUE) %>% mutate(next_start = lead(start)) 
  # Result: 
  # (a*)|exon3|(a)----<------(e*)|exon 2|(e)----   ---> dataframe of 3 columns for this isoform: start = a*, end = a, next_start = e* 
  
  # For Method 2, A5' or A3' have to differ between the start coordinates of the exon, but the same for end coordinates the previous exon from isoform, and same end coordinates of the same exon
  # i.e Isoform 3 and Isoform 1, whereby abs(c-b) == 0 but abs(g*-f*) != 0
  # To calculate difference between start coordinates (c-b) => lag the start coordinates
  # To calculate difference between end coordinates from the previous exon from isoform
  input_dat <- input_dat %>% group_by(isoform) %>% arrange(start, .by_group = TRUE) %>% mutate(previous_end = lag(end))
  # Result: 
  # (c*)|exon1|(c)---->-----------(f*)|exon 2|(f)----  ---> dataframe of 3 columns for this isoform: start = f*, end = f, previous_end = c* 
  
  # While A5' and A3' cannot occur from final or first exon, by making sure that both conditions match (diff at end, and diff at next_start), any cases arising from difference in exonic length due to first or last exon would be filtered out 
  #1        (b*)|exon20|(b)---->------(e*)|exon 21 (final)|(e)
  #2        (c*)|exon20|(c)---->------(f*)|exon 21 (final         |(f)
  #In example above, Isoform 1 and Isoform 2 may appear to have 5' if only using diff_at_end (f - e), however next_start would be NA as this is the final exon; therefore this would be filtered out 
  
  
  Diff_coordinates <- input_dat %>% group_by(associated_gene) %>% arrange(end, .by_group = TRUE) %>% 
    mutate(Diff_end = abs(end - lag(end))) %>%
    mutate(Diff_start = abs(start - lag(start))) %>%
    mutate(Diff_next_start = abs(next_start - lag(next_start))) %>% 
    mutate(Diff_previous_end = abs(previous_end - lag(previous_end))) 
  
  # THRESHOLD!
  # standard exon length is 300
  #Method1_wobble_hist <- Diff_coordinates %>% filter(Diff_end < 1000, Diff_end > 0) %>% ggplot(., aes(x = Diff_end)) + geom_histogram(bins = 50) + theme_bw() + 
  #labs(x = "Difference in end coordinates per exon across isoforms")
  #Method2_wobble_hist <- Diff_coordinates %>% filter(Diff_start < 1000, Diff_start > 0) %>% ggplot(., aes(x = Diff_start)) + geom_histogram(bins = 50) + theme_bw() + 
  #labs(x = "Difference in end coordinates per exon across isoforms")
  # Method1_wobble_hist and Method2_wobble_hist very similar, with hist falling after 250bp, referring to the distance of the next exon rather than A5'/A3' (take 200bp)
  
  Method1 <- Diff_coordinates %>% 
    # Not A3'/A5' if Diff_end is NA (i.e first entry for each isoform), no difference in end, Diff_end > threshold specified (so next exon along), or if the start is not the same, or wobble of 2bp 
    mutate(Diff_class_end = ifelse(is.na(Diff_end) | Diff_end == 0 | Diff_end > threshold | Diff_start != 0 | Diff_end < 3, "No","Yes")) %>%
    mutate(Diff_class_start = ifelse(Diff_next_start == 0, "Yes","No")) %>% 
    # Remove Diff_previous_end as this refers to detected lengths due to first exon (alternative first)
    filter(Diff_class_end == "Yes" & Diff_class_start == "Yes" & !is.na(Diff_previous_end)) %>% 
    mutate(Event = ifelse(strand == "+", "A5", "A3"))
  
  Method2 <- Diff_coordinates %>%
    mutate(Diff_class_start = ifelse(is.na(Diff_start) | Diff_start == 0 | Diff_start > threshold | Diff_end != 0 | Diff_start < 3 , "No","Yes")) %>%
    mutate(Diff_class_end = ifelse(Diff_previous_end == 0, "Yes","No")) %>% 
    filter(Diff_class_end == "Yes" & Diff_class_start == "Yes" & !is.na(Diff_next_start)) %>% 
    mutate(Event = ifelse(strand == "+", "A3", "A5"))
  
  Method1 <<- Method1 
  Method2 <<- Method2
  A3_A5 <- rbind(Method1, Method2)
  return(A3_A5)
}


