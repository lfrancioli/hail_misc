#Load data

indels = read.delim("introns_counts.txt.bgz", stringsAsFactors = F)

ggplot(indels,aes(nIndels)) + geom_density()

indels %<>%
  mutate(
    pass = pass == "true",
    lcr = lcr == "true",
    segdup = segdup == "true",
    wasSplit = wasSplit == "true",
    intron_length = end-start,
    indel_per_bp = nIndels / intron_length,
    indel_middle_per_bp = (nIndels - nIndels_200bp_start - nIndels_200bp_end) / (intron_length - 400),
    indels_ends = nIndels_200bp_start + nIndels_200bp_end
    )

all_pass_indels = indels %>%
  filter(pass) %>%
  group_by(intron_id) %>%
  summarise(
    metric_entropy = first(metric_entropy),
    start_200bp_entropy = first(start_200bp_entropy),
    end_200bp_entropy = first(end_200bp_entropy),
    ends_entropy = first(ends_entropy),
    middle_entropy = first(middle_entropy),
    nIndels = sum(nIndels),
    nIndels_200bp_start = sum(nIndels_200bp_start),
    nIndels_200bp_end = sum(nIndels_200bp_end),
    intron_length = first(intron_length),
    indel_per_bp = nIndels / intron_length,
    indel_middle_per_bp = (nIndels - nIndels_200bp_start - nIndels_200bp_end) / (intron_length - 400),
    indels_ends = nIndels_200bp_start + nIndels_200bp_end
  )

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
           data=filtered_indels))

#Look at all pass indels
summary(lm(indel_per_bp ~ intron_length + metric_entropy + intron_length*metric_entropy,
           data=all_pass_indels))

#Limit to introns < 120bp
summary(lm(indel_per_bp ~ small_intron + metric_entropy + intron_length*metric_entropy,
           data=all_pass_indels %>% mutate(
             small_intron = intron_length < 100
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
  