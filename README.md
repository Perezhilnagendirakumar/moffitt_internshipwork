# "Deciphering Cancer Heterogeneity: Computational Biology Insights into immuno microenvironment differences in Head and Neck Squamous Cell Carcinoma (HNSCC) ".

The introduction of single-cell RNA sequencing (scRNA-seq) technology has transformed biomedical research, allowing the detailed study of individual cells. The increasing data throughput and enhanced efficiency have resulted in scRNA-seq datasets containing transcription profiles of over a million cells, accompanied by a substantial increase in available analysis tools. Despite these advancements, fully dissecting cellular heterogeneity within a large cell population remains a computational challenge. This repo presents steps needed to make sense of single-cell RNA sequencing (scRNA) data. I used a scRNA dataset from the article Immune landscape of viral- and carcinogen-drived head and neck cancer published in Nature Communications. 

### ScRNA analysis using Seurat  
Single-cell RNA sequencing (scRNA-seq) datasets were downloaded from the Gene Expression Omnibus through accession number GEO: GSE139324.The study included measurements from 131,224 individual immune cells in the TME from patients with HPV− (n = 18) and HPV + HNSCC (n = 8) who were immunotherapy treatment-naïve. Seurat V5 was used for preprocessing , Integration , batch effect reduction and unsupervised clustering . 
 
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/dc1700b4-aa9d-4b8b-9a6f-9ad4a0e2d472)![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/d2d130be-7803-4937-a6d6-be478978b136)![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/87089ce9-39f6-4f0d-905e-e86a63223566)![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/f5f02d55-458b-4225-b20e-6d1783897243)


###  Annotation of cell type using marker genes  and ultrafast cell type clssification tool "SC type"

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/02385d16-3e35-4eaf-a593-7bf6589d33a2)

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/3061e817-8115-4b92-93a7-cd04507b7543)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/e69d99cd-5416-41ba-9725-b54fbe4e85df)![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/b3146aa8-c6aa-4131-8ba3-5a293bc6367e)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/d78332c4-ddad-4b83-ab77-1ae4570ff1e6)

###  Exploring HPV + enriched Cell type by subclustering and reanalysis using Seurat 

According to the Rational Barplot mentioned above two group of cell types like B cells and Naive T cells are highly enriched in HPV + tumors . Among which B cell population were higher when compared to Naive T cells  . Performed Reclustering and analysis of B cells using seurat and annotated their cell type based on DEG 


![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/d0ddb107-32d5-4f13-b0b7-5a479a8c677d)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/fad8042e-5802-4c94-85b6-8063f0b1fc0b)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/b6660e45-9e9f-40a2-8d71-0626ba67e202)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/2f3425ec-d28f-420c-9c5d-d6580e84fdce)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/fe13a1ad-9d7f-424c-bee4-b61d1002ddfb)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/0e5598c1-2074-49e6-8680-1eeabec353ca)![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/85d8b73e-e851-4e27-91c9-7f5b3747c44b)



###  Enrichment analysis using Clusterprofiler 

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/dba4e504-78a5-41e2-a68c-cbcd06f0dd53)

**GO_ BP**

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/ae645747-e68d-4ac5-8d19-57fb90df308a)

**KEGG-BP**

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/37e61313-9e41-4477-af67-0c928913d9c5)

###  Pseudotimetrajectory analysis using monocle3
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/abb98c42-de37-407d-8c68-6873c0f993b2)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/120ae284-7e98-4225-a6ff-39f7b9c3e444)


### Rare cell /minor cell population  Identification using FIRE, Giniclust2 and Gapclust 

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/e13d136b-5234-496e-b2a2-69f0cdf2a697)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/c0703de3-ff5f-4982-b576-e3c96f998011)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/571177c4-6ecd-4c0f-8ccc-c522e52aff8d)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/03f33114-33d9-4017-8293-beba9e3bd220)

### WGCNA (Weighted Gene Coexpression network analysis) 
Clustering dendrogram of genes. Gene clustering tree (dendrogram) obtained by hierarchical clustering of adjacency-based dissimilarity. The colored row below the dendrogram indicates module membership identified by the dynamic tree cut method, together with assigned merged module colors and the original module colors.

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/03405cfc-791e-4983-8427-19df3bdbc9c8) 
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/cf92ae9d-e9fa-4616-ae18-9e3bf4451604)


### Identification of the modules of interest that best represent the gene expression in rare cell or minor cell populations identified by tools such as FIRE, Giniclust2, and GapClust.

![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/08a7014d-b4a0-422b-86be-3c2fe3de65a3)
![image](https://github.com/Perezhilnagendirakumar/moffitt_internshipwork/assets/97453603/2c9be4cc-deee-46aa-8246-bfdce6def0bb)


###  Conclusion 

The gene expression profiles of minor cell populations were closely associated with Type I interferon cellular responses, as evidenced by Gene Ontology Biological Process (GO_BP) enrichment analysis. Type I interferons (IFNs) are pivotal cytokines in antiviral defense, orchestrating both innate and adaptive immune responses. IFNs amplify B lymphocyte responses to viral challenges, fostering the production of cytotoxic and neutralizing antibodies. Their multifaceted impact on B-cell function encompasses receptor engagement, Toll-like receptor expression, cell migration, antigen presentation, cytokine responsiveness, cytokine production, survival, differentiation, and class-switch recombination.Remarkably, the Type I interferon response was significantly enriched in HPV+ HNSCC populations. This Type I interferon response was prominently observed during the transition from activated B cells(ABC)to antibody-secreting B cells (ASC), peaking as pseudotime increased, and gradually declining post-transition. This dynamic regulation underscores the role of Type I interferons in the differentiation process of B cells, further emphasizing their relevance in the immune landscape of HPV+ HNSCC and potential therapeutic avenues.


























































