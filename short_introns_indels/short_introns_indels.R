#Load data

ends_length = 40
introns_gencode = read.delim(sprintf("introns_counts.%dbp.txt.bgz", ends_length), stringsAsFactors = F)
introns = read.delim(sprintf("introns_count.%dbp.refseq.txt", ends_length), stringsAsFactors = F)


ggplot(indels,aes(nIndels)) + geom_density()

introns %<>%
  replace_na(list(n = 0, n_start = 0, n_end =0)) %>%
  mutate(
    pass = pass == "true",
    lcr = lcr == "true",
    segdup = segdup == "true",
    wasSplit = wasSplit == "true",
    deletion = deletion == "true",
    intron_length = end-start,
    n_per_bp = n / intron_length,
    n_middle_per_bp = (n - n_start - n_end) / (intron_length - 2*ends_length),
    n_ends = n_start + n_end
    ) %>%
  replace_na(list(pass = T, lcr = F, segdup = F, wasSplit = F, deletion = F))

introns_pass_indels = introns %>%
  filter(pass) %>%
  group_by(intron_id) %>%
  summarise(
    metric_entropy = first(metric_entropy),
    start_entropy = first(start_entropy),
    end_entropy = first(end_entropy),
    ends_entropy = first(ends_entropy),
    middle_entropy = first(middle_entropy),
    n_indels = sum(n),
    n_deletions = sum(n[deletion]),
    n_insertions = n_indels - n_deletions,
    in_del_ratio = n_insertions / n_deletions,
    n_indels_start = sum(n_start),
    n_deletions_start = sum(n_start[deletion]),
    n_insertions_start = n_indels_start - n_deletions_start,
    n_indels_end = sum(n_end),
    n_deletions_end = sum(n_end[deletion]),
    n_insertions_end = n_indels_end - n_deletions_end,
    n_indels_ends = n_indels_start + n_indels_end,
    n_deletions_ends = n_deletions_start + n_deletions_end,
    n_insertions_ends = n_insertions_start + n_insertions_end,
    intron_length = first(intron_length),
    indel_per_bp = n_indels / intron_length,
    deletion_per_bp = n_deletions / intron_length,
    insertion_per_bp = n_insertions / intron_length,
    indel_middle_per_bp = (n_indels - n_indels_start - n_indels_end) / (intron_length - 400),
    indels_ends = n_indels_start + n_indels_end,
    small_intron = intron_length < 2*ends_length
  ) %>% 
  rowwise() %>%
  mutate(
    ends_length = min(ends_length, intron_length),
    insertions_start_per_bp = n_insertions_start / min(ends_length/2, intron_length),
    insertions_end_per_bp = n_insertions_end / min(ends_length/2, intron_length),
    insertions_ends_per_bp = n_insertions_ends / min(ends_length, intron_length),
    deletions_start_per_bp = n_deletions_start / min(ends_length/2, intron_length),
    deletions_end_per_bp = n_deletions_end / min(ends_length/2, intron_length),
    deletions_ends_per_bp = n_deletions_ends / min(ends_length, intron_length),
    indels_start_per_bp = n_indels_start / min(ends_length/2, intron_length),
    indels_end_per_bp = n_indels_end / min(ends_length/2, intron_length),
    indels_ends_per_bp = n_indels_ends / min(ends_length, intron_length)
  ) %>%
  ungroup() %>% 
  filter(intron_length > 10)

##Plot entropy distribution for first 100bp of introns
introns_pass_indels %>%
  ggplot(aes(start_entropy, col=small_intron)) +
  geom_density()

#Test if significant difference  
wilcox.test(introns_pass_indels %>% filter(small_intron) %$% start_entropy,
            introns_pass_indels %>% filter(!small_intron) %$% start_entropy)
ks.test(introns_pass_indels %>% filter(small_intron) %$% start_entropy,
            introns_pass_indels %>% filter(!small_intron) %$% start_entropy)

#Look at tabled numbers
introns_pass_indels_summary = introns_pass_indels %>%
  group_by(small_intron) %>%
  summarise(n = n(),
            n_ins = sum(insertions_ends_per_bp * (2*ends_length)),
            n_del = sum(deletions_ends_per_bp * (2*ends_length)),
            n_has_ins = sum(n_insertions_ends > 0),
            n_has_del = sum(n_deletions_ends > 0),
            tot_ins_per_bp = sum(n_insertions_ends) / sum(ends_length),
            tot_del_per_bp = sum(n_deletions_ends) / sum(ends_length),
            prop_has_ins = n_has_ins / n(),
            prop_has_del = n_has_del / n()
            )

#Test for number of del vs ins
fisher.test(
matrix(c(introns_pass_indels_summary %$% n_ins,
         introns_pass_indels_summary %$% n_del),
       nrow = 2,
       dimnames = list(c("long intron","short intron"),
                      c("n_ins","n_del"))
)
)

#Test for presence / absence of dels
fisher.test(
  matrix(c(introns_pass_indels_summary %$% n - introns_pass_indels_summary %$% n_has_del,
           introns_pass_indels_summary %$% n_has_del),
         nrow = 2,
         dimnames = list(c("long intron","short intron"),
                         c("n_no_del","n_del"))
  )
)

 ggplot(introns_pass_indels, aes(intron_length, in_del_ratio)) + 
  geom_point() +
  scale_x_log10()

ggplot(indels,aes(intron_length)) + geom_histogram() + scale_x_log10()

ggplot(indels,aes(indel_per_bp)) + geom_density()
ggplot(indels,aes(metric_entropy)) + geom_density()

#Filtered to best indels in introns > 500bp
filtered_indels = indels %>%
  filter(pass & 
         !lcr &
         !segdup &
         !wasSplit)


#Filter to most reliable indels
summary(lm(indel_per_bp ~ intron_length + metric_entropy + intron_length*metric_entropy,
           data=introns_pass_indels))

#Look at all pass indels
summary(lm(indel_per_bp ~ intron_length + metric_entropy + intron_length*metric_entropy,
           data=all_pass_indels))


ggplot(introns_pass_indels, aes(in_del_ratio, col = intron_length < 60)) + geom_density()
ggplot(introns_pass_indels, aes(insertion_per_bp, col = intron_length < 60)) + geom_density()



introns_pass_indels %>% 
  mutate(small_intron = intron_length < 80) %>% 
  group_by(small_intron) %>% 
  summarise(n=n(), n_no_del = sum(n_deletions ==0 ), n_dels = sum(n_deletions), n_no_ins = sum(n_insertions ==0 ), n_ins = sum(n_insertions))


introns_pass_indels %>% 
  rowwise() %>%
  mutate(
    small_intron = intron_length < 80,
    deletion_per_bp_start = n_deletions_start / min(100, intron_length)
  ) %>%
  ggplot(aes(deletion_per_bp_start, col=small_intron)) + geom_density()

#Look at small introns
summary(lm(deletion_per_bp_start ~ small_intron + metric_entropy,
           data=introns_pass_indels %>% 
             rowwise() %>%
             mutate(
             small_intron = intron_length < 80,
             deletion_per_bp_start = n_deletions_start / min(100, intron_length)
           )))

summary(lm(insertion_per_bp_start ~ small_intron + metric_entropy,
           data=introns_pass_indels %>% 
             rowwise() %>%
             mutate(
               small_intron = intron_length < 80,
               insertion_per_bp_start = n_insertions_start / min(100, intron_length)
             )))

#look at best indels -- middle part only
summary(lm(indel_middle_per_bp ~ intron_length + middle_entropy,
           data = indels %>% filter(pass & intron_length > 500)))

#look at all pass indels -- middle part only
summary(lm(indel_middle_per_bp ~ intron_length + middle_entropy,
        data = all_pass_indels %>% filter(intron_length > 500)))

#Look at ends of introns only
summary(lm(indels_ends ~ intron_length + ends_entropy),
        data )

#Best quality indels only
summary(lm(indel_middle_per_bp ~ intron_length + metric_entropy,
           data = filtered_indels %>% 
             filter(intron_length > 500)))

#Plot this
filtered_indels %>%
  ggplot(aes(indel_middle_per_bp, col = intron_length > 1000)) +
  geom_density()

filtered_indels %>%
  ggplot(aes(start_200bp_entropy + end_200bp_entropy, col = intron_length > 1000)) +
  geom_density()

filtered_indels %>%
  ggplot(aes(metric_entropy, col = intron_length > 1000)) +
  geom_density()

filtered_indels %>%
  ggplot(aes(nIndels_200bp_start, col = intron_length < 1000)) +
  geom_density()

summary(lm(nIndels_200bp_start ~ intron_length + start_200bp_entropy,
           data = filtered_indels))
  