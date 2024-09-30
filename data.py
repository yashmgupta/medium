# data.py

research_summaries = [
    # 1. Biological Sciences

    # Molecular Biology
    {
        "id": 1,
        "category": "Biological Sciences",
        "title": "Molecular Biology - DNA, RNA, and Protein Structure and Function",
        "content": """
**DNA, RNA, and Protein Structure and Function**

**DNA (Deoxyribonucleic Acid):**
- **Structure:** Double helix composed of nucleotides (adenine, thymine, cytosine, guanine).
- **Function:** Stores genetic information crucial for the development, functioning, and reproduction of all living organisms.

**RNA (Ribonucleic Acid):**
- **Structure:** Typically single-stranded, composed of nucleotides (adenine, uracil, cytosine, guanine).
- **Function:** Plays roles in coding, decoding, regulation, and expression of genes. Types include mRNA, tRNA, and rRNA.

**Proteins:**
- **Structure:** Polypeptide chains folded into specific three-dimensional shapes, consisting of amino acids.
- **Function:** Serve as enzymes, structural components, signaling molecules, and perform a myriad of cellular functions.

**Central Dogma of Molecular Biology:**
- **Flow of Genetic Information:** DNA → RNA → Protein.
- **Processes Involved:**
  - **Transcription:** Synthesis of RNA from DNA template.
  - **Translation:** Synthesis of proteins from mRNA template.
""",
        "tags": ["Molecular Biology", "DNA", "RNA", "Protein Structure", "Genetics"],
    },
    {
        "id": 2,
        "category": "Biological Sciences",
        "title": "Molecular Biology - Central Dogma of Molecular Biology",
        "content": """
**Central Dogma of Molecular Biology**

**Definition:**
- Describes the flow of genetic information within a biological system.
- **Primary Pathway:** DNA → RNA → Protein.

**Processes:**
1. **Replication:**
   - **Description:** DNA makes a copy of itself during cell division.
   - **Enzymes Involved:** DNA polymerase, helicase.

2. **Transcription:**
   - **Description:** DNA is transcribed into messenger RNA (mRNA).
   - **Enzymes Involved:** RNA polymerase.
   - **Steps:**
     - Initiation at promoter regions.
     - Elongation of the RNA strand.
     - Termination upon reaching termination signals.

3. **Translation:**
   - **Description:** mRNA is translated into a specific protein sequence.
   - **Components Involved:** Ribosomes, transfer RNA (tRNA), amino acids.
   - **Steps:**
     - Initiation: Ribosome assembles around the start codon.
     - Elongation: tRNAs bring amino acids to the ribosome.
     - Termination: Ribosome reaches a stop codon, releasing the protein.
""",
        "tags": ["Central Dogma", "Transcription", "Translation", "Gene Expression"],
    },
    {
        "id": 3,
        "category": "Biological Sciences",
        "title": "Molecular Biology - Transcription and Translation",
        "content": """
**Transcription and Translation**

**Transcription:**
- **Purpose:** Convert DNA sequences into RNA.
- **Location:** Nucleus (in eukaryotes).
- **Steps:**
  1. **Initiation:** RNA polymerase binds to promoter region of DNA.
  2. **Elongation:** RNA strand is synthesized complementary to DNA template.
  3. **Termination:** RNA synthesis stops at termination signal; RNA molecule is released.

**Translation:**
- **Purpose:** Synthesize proteins based on mRNA sequence.
- **Location:** Cytoplasm, on ribosomes.
- **Steps:**
  1. **Initiation:** Ribosome assembles around the start codon (AUG) on mRNA.
  2. **Elongation:** tRNA molecules bring amino acids corresponding to codons; ribosome links amino acids into a polypeptide chain.
  3. **Termination:** Ribosome encounters a stop codon (UAA, UAG, UGA); protein synthesis ends, and the ribosome disassembles.

**Key Components:**
- **mRNA (Messenger RNA):** Carries genetic information from DNA to ribosome.
- **tRNA (Transfer RNA):** Brings amino acids to the ribosome during protein synthesis.
- **Ribosomes:** Molecular machines that facilitate translation by reading mRNA and assembling proteins.
""",
        "tags": ["Transcription", "Translation", "Protein Synthesis", "mRNA", "tRNA", "Ribosomes"],
    },

    # Genetics
    {
        "id": 4,
        "category": "Biological Sciences",
        "title": "Genetics - Mendelian Genetics",
        "content": """
**Mendelian Genetics**

**Gregor Mendel:**
- Known as the "Father of Genetics."
- Conducted experiments on pea plants to understand inheritance patterns.

**Key Principles:**
1. **Law of Segregation:**
   - Each organism carries two alleles for each trait.
   - Alleles segregate during gamete formation, ensuring each gamete carries only one allele.

2. **Law of Independent Assortment:**
   - Alleles of different genes assort independently during gamete formation.
   - Applies to genes located on different chromosomes or far apart on the same chromosome.

**Genotypic and Phenotypic Ratios:**
- **Monohybrid Cross:** Involves one trait.
  - **Genotypic Ratio:** 1 homozygous dominant : 2 heterozygous : 1 homozygous recessive.
  - **Phenotypic Ratio:** 3 dominant phenotype : 1 recessive phenotype.

- **Dihybrid Cross:** Involves two traits.
  - **Genotypic Ratio:** 9:3:3:1.
  - **Phenotypic Ratio:** 9 dominant for both traits : 3 dominant for one and recessive for the other : 3 recessive for one and dominant for the other : 1 recessive for both traits.

**Punnett Squares:**
- Tool used to predict the probability of offspring genotypes and phenotypes.
""",
        "tags": ["Genetics", "Mendelian Genetics", "Inheritance", "Punnett Square"],
    },
    {
        "id": 5,
        "category": "Biological Sciences",
        "title": "Genetics - Population Genetics",
        "content": """
**Population Genetics**

**Definition:**
- Study of genetic variation within populations and involves the examination of allele frequency changes under influences such as mutation, natural selection, genetic drift, gene flow, and mating patterns.

**Key Concepts:**
1. **Hardy-Weinberg Equilibrium:**
   - Principle stating that allele and genotype frequencies in a population will remain constant from generation to generation in the absence of evolutionary influences.
   - **Conditions:** No mutation, random mating, no gene flow, infinite population size, and no selection.

2. **Genetic Drift:**
   - Random fluctuations in allele frequencies, especially in small populations.
   - Can lead to the loss or fixation of alleles over time.

3. **Gene Flow:**
   - Transfer of alleles or genes from one population to another.
   - Increases genetic diversity within a population.

4. **Mutation:**
   - Random changes in the DNA sequence.
   - Introduces new genetic variation into a population.

5. **Natural Selection:**
   - Differential survival and reproduction of individuals due to differences in phenotype.
   - Leads to adaptation of populations to their environments.

**Genetic Variation:**
- **Sources:** Mutation, sexual reproduction, gene flow.
- **Importance:** Provides the raw material for evolution and adaptation.

**Applications:**
- Understanding evolutionary processes.
- Conservation genetics.
- Medical genetics and understanding genetic diseases.
""",
        "tags": ["Population Genetics", "Allele Frequency", "Genetic Drift", "Gene Flow", "Hardy-Weinberg"],
    },
    {
        "id": 6,
        "category": "Biological Sciences",
        "title": "Genetics - Quantitative Genetics",
        "content": """
**Quantitative Genetics**

**Definition:**
- Branch of genetics that deals with phenotypes that vary continuously (e.g., height, weight) and are typically influenced by multiple genes and environmental factors.

**Key Concepts:**
1. **Polygenic Traits:**
   - Traits controlled by multiple genes, often exhibiting a range of phenotypes.
   - Examples: Height, skin color, intelligence.

2. **Heritability:**
   - Measure of how much of the variation in a trait is due to genetic factors versus environmental factors.
   - **Formula:** \( h^2 = \frac{V_G}{V_P} \)
     - \( V_G \): Genetic variance.
     - \( V_P \): Phenotypic variance.

3. **Additive Genetic Variance:**
   - Part of genetic variance attributed to the additive effect of different alleles.
   - Predicts response to selection.

4. **Non-Additive Genetic Variance:**
   - Includes dominance and epistatic variance.
   - Less predictable in response to selection.

5. **Selection:**
   - **Directional Selection:** Favors one extreme phenotype.
   - **Stabilizing Selection:** Favors intermediate phenotypes.
   - **Disruptive Selection:** Favors both extreme phenotypes.

**Applications:**
- Breeding programs in agriculture.
- Understanding the genetic basis of complex diseases.
- Evolutionary biology studies.
""",
        "tags": ["Quantitative Genetics", "Polygenic Traits", "Heritability", "Genetic Variance"],
    },

    # Genomics
    {
        "id": 7,
        "category": "Biological Sciences",
        "title": "Genomics - Genome Organization and Evolution",
        "content": """
**Genome Organization and Evolution**

**Genome Organization:**
- **Prokaryotic Genomes:**
  - Typically circular DNA molecules.
  - Contain essential genes for survival and reproduction.
  - Often include plasmids (small, circular DNA molecules) that carry non-essential genes.

- **Eukaryotic Genomes:**
  - Linear chromosomes housed within a nucleus.
  - Contain a large number of non-coding regions (introns, regulatory sequences).
  - Organized into chromatin (DNA-protein complex).

**Genome Size and Complexity:**
- Not directly correlated with organismal complexity.
- **C-value Paradox:** Larger genomes do not necessarily mean more complex organisms.
- Variation in genome size due to differences in non-coding DNA, repetitive elements, and gene duplication events.

**Genome Evolution:**
1. **Gene Duplication:**
   - Duplication of genes leading to genetic redundancy.
   - Provides raw material for the evolution of new functions.

2. **Horizontal Gene Transfer:**
   - Transfer of genes between different species.
   - Common in prokaryotes, contributing to genetic diversity.

3. **Mutation and Selection:**
   - Accumulation of mutations introduces genetic variation.
   - Natural selection acts on this variation, shaping genomes over time.

4. **Transposable Elements:**
   - DNA sequences that can change their position within the genome.
   - Contribute to genetic diversity and genome plasticity.

5. **Genome Rearrangements:**
   - Large-scale structural changes like inversions, translocations, deletions.
   - Can have significant impacts on gene function and regulation.

**Comparative Genomics:**
- Study of the similarities and differences in the genomes of different species.
- Helps in understanding evolutionary relationships and functional genomics.
""",
        "tags": ["Genomics", "Genome Organization", "Genome Evolution", "Comparative Genomics"],
    },
    {
        "id": 8,
        "category": "Biological Sciences",
        "title": "Genomics - Comparative Genomics",
        "content": """
**Comparative Genomics**

**Definition:**
- Branch of genomics that involves comparing the genomes of different species to identify similarities and differences.
- Helps in understanding the structure, function, and evolutionary relationships of genomes.

**Key Areas:**
1. **Orthologous Genes:**
   - Genes in different species that originated from a common ancestral gene through speciation.
   - Typically retain the same function.

2. **Paralogous Genes:**
   - Genes related by duplication within a genome.
   - May evolve new functions.

3. **Synteny:**
   - Conservation of blocks of genes across different species.
   - Indicates evolutionary conservation and functional importance.

4. **Genome Synteny Analysis:**
   - Identifies conserved gene order and orientation.
   - Useful for genome assembly and annotation.

**Applications:**
- **Functional Annotation:**
  - Predicting gene function based on conserved sequences across species.

- **Evolutionary Biology:**
  - Tracing the evolutionary history of genes and genomes.
  - Understanding mechanisms of speciation and adaptation.

- **Identifying Conserved Elements:**
  - Regulatory regions, non-coding RNAs, and essential genes conserved across species.

- **Comparative Mapping:**
  - Aligning genetic maps of different species to identify conserved regions.

**Tools and Databases:**
- **Ensembl:**
  - Genome browser providing access to comparative genomics data.

- **UCSC Genome Browser:**
  - Offers tools for visualizing synteny and comparative data.

- **OrthoDB:**
  - Database of orthologous genes across multiple species.

- **MAUVE:**
  - Software for multiple genome alignment.
""",
        "tags": ["Comparative Genomics", "Orthologous Genes", "Paralogous Genes", "Synteny"],
    },
    {
        "id": 9,
        "category": "Biological Sciences",
        "title": "Genomics - Functional Genomics",
        "content": """
**Functional Genomics**

**Definition:**
- Area of genomics that aims to describe gene functions and interactions.
- Integrates data from various omics technologies to understand gene expression, regulation, and interaction networks.

**Key Areas:**
1. **Gene Expression Analysis:**
   - Studying when and where genes are active.
   - Techniques include RNA sequencing (RNA-seq) and microarrays.

2. **Regulatory Genomics:**
   - Understanding how gene expression is controlled.
   - Focuses on regulatory elements like promoters, enhancers, and transcription factors.

3. **Interaction Networks:**
   - Mapping interactions between proteins, genes, and other molecules.
   - Utilizes tools like yeast two-hybrid screens and co-immunoprecipitation.

4. **Epigenomics:**
   - Study of heritable changes in gene expression that do not involve changes to the DNA sequence.
   - Includes DNA methylation, histone modification, and chromatin remodeling.

5. **Comparative Functional Genomics:**
   - Comparing gene functions across different species to identify conserved and divergent functions.

**Techniques and Technologies:**
- **Chromatin Immunoprecipitation Sequencing (ChIP-seq):**
  - Identifies binding sites of DNA-associated proteins.
  
- **RNA Interference (RNAi) and CRISPR-Cas9:**
  - Tools for gene knockdown and knockout to study gene function.

- **Mass Spectrometry:**
  - Identifies and quantifies proteins to understand the proteome.

**Applications:**
- **Disease Research:**
  - Identifying genes involved in diseases and understanding their mechanisms.

- **Drug Discovery:**
  - Target identification and validation through understanding gene functions.

- **Agricultural Genomics:**
  - Enhancing crop traits by understanding and manipulating gene functions.

- **Synthetic Biology:**
  - Designing and constructing new biological parts and systems based on functional genomics insights.
""",
        "tags": ["Functional Genomics", "Gene Expression", "Regulatory Genomics", "Epigenomics"],
    },

    # Proteomics
    {
        "id": 10,
        "category": "Biological Sciences",
        "title": "Proteomics - Protein Structure and Function",
        "content": """
**Protein Structure and Function**

**Protein Structure:**
1. **Primary Structure:**
   - Linear sequence of amino acids in a polypeptide chain.
   - Determined by the gene encoding the protein.

2. **Secondary Structure:**
   - Localized folding patterns within a protein, stabilized by hydrogen bonds.
   - Common structures include alpha-helices and beta-sheets.

3. **Tertiary Structure:**
   - Overall three-dimensional shape of a single protein molecule.
   - Formed by interactions between secondary structures and side chains of amino acids.

4. **Quaternary Structure:**
   - Arrangement of multiple protein subunits in a multi-subunit complex.
   - Example: Hemoglobin, which consists of four subunits.

**Protein Function:**
- **Enzymes:** Catalyze biochemical reactions.
- **Structural Proteins:** Provide support and shape to cells and tissues.
- **Transport Proteins:** Carry molecules across cell membranes or within the bloodstream.
- **Signaling Proteins:** Involved in cell communication and signal transduction.
- **Immune Proteins:** Protect against pathogens (e.g., antibodies).

**Structure-Function Relationship:**
- The specific three-dimensional structure of a protein determines its function.
- **Active Sites:** Regions on enzymes where substrates bind and reactions occur.
- **Binding Sites:** Areas where proteins interact with other molecules or proteins.

**Protein Domains:**
- Distinct functional and/or structural units within a protein.
- Can fold independently and contribute to the protein's overall function.

**Post-Translational Modifications (PTMs):**
- Chemical modifications that occur after protein synthesis.
- Examples include phosphorylation, glycosylation, ubiquitination.
- PTMs can alter protein activity, localization, stability, and interactions.
""",
        "tags": ["Proteomics", "Protein Structure", "Protein Function", "Enzymes", "PTMs"],
    },
    {
        "id": 11,
        "category": "Biological Sciences",
        "title": "Proteomics - Post-Translational Modifications",
        "content": """
**Post-Translational Modifications (PTMs)**

**Definition:**
- Chemical modifications that occur to proteins after their synthesis (translation).
- Regulate protein function, activity, localization, and interactions.

**Common Types of PTMs:**
1. **Phosphorylation:**
   - Addition of a phosphate group, typically to serine, threonine, or tyrosine residues.
   - Regulates enzyme activity, signal transduction, and protein interactions.

2. **Glycosylation:**
   - Attachment of carbohydrate chains to asparagine (N-linked) or serine/threonine (O-linked) residues.
   - Affects protein folding, stability, and cell-cell recognition.

3. **Ubiquitination:**
   - Addition of ubiquitin proteins to lysine residues.
   - Marks proteins for degradation via the proteasome.

4. **Methylation:**
   - Addition of methyl groups to lysine or arginine residues.
   - Involved in regulating gene expression and protein-protein interactions.

5. **Acetylation:**
   - Addition of acetyl groups to lysine residues.
   - Influences gene expression by modifying histones and affects protein stability.

6. **Lipidation:**
   - Attachment of lipid molecules to proteins.
   - Facilitates membrane association and protein localization.

7. **Sumoylation:**
   - Addition of Small Ubiquitin-like Modifier (SUMO) proteins.
   - Modulates protein function, localization, and interactions.

**Functions and Implications of PTMs:**
- **Regulation of Enzyme Activity:**
  - Phosphorylation can activate or inhibit enzymes.
  
- **Signal Transduction:**
  - PTMs play critical roles in transmitting signals within and between cells.
  
- **Protein Stability and Degradation:**
  - Ubiquitination targets proteins for degradation, controlling protein levels.
  
- **Subcellular Localization:**
  - PTMs can direct proteins to specific cellular compartments.
  
- **Protein-Protein Interactions:**
  - Modifications can enhance or inhibit interactions between proteins.

**Detection and Analysis of PTMs:**
- **Mass Spectrometry (MS):**
  - Widely used for identifying and quantifying PTMs.
  
- **Western Blotting:**
  - Uses antibodies specific to certain PTMs.
  
- **Enzyme-Linked Immunosorbent Assay (ELISA):**
  - Quantifies PTMs using specific antibodies.

**Biological Significance:**
- PTMs are essential for the dynamic regulation of cellular processes.
- Dysregulation of PTMs is associated with various diseases, including cancer, neurodegenerative disorders, and diabetes.
""",
        "tags": ["Proteomics", "Post-Translational Modifications", "Phosphorylation", "Ubiquitination"],
    },
    {
        "id": 12,
        "category": "Biological Sciences",
        "title": "Proteomics - Protein-Protein Interactions",
        "content": """
**Protein-Protein Interactions (PPIs)**

**Definition:**
- Physical contacts established between two or more protein molecules as a result of biochemical events and/or electrostatic forces.
- PPIs are fundamental to most biological processes.

**Types of PPIs:**
1. **Transient Interactions:**
   - Short-lived and often involved in signaling pathways.
   - Example: Enzyme-substrate complexes.

2. **Permanent Interactions:**
   - Stable associations forming multi-protein complexes.
   - Example: Hemoglobin, which consists of four subunits.

3. **Obligate Interactions:**
   - Proteins require each other to maintain their structure.
   - Example: α and β subunits of hemoglobin.

4. **Non-Obligate Interactions:**
   - Proteins can exist independently but interact under certain conditions.
   - Example: Transcription factors binding to DNA.

**Methods to Detect PPIs:**
1. **Yeast Two-Hybrid Screening:**
   - Detects binary interactions between proteins.
   - Based on reconstitution of a transcription factor in yeast.

2. **Co-Immunoprecipitation (Co-IP):**
   - Uses antibodies to pull down a target protein and its binding partners from a cell lysate.

3. **Affinity Purification Mass Spectrometry (AP-MS):**
   - Combines affinity purification with mass spectrometry to identify interacting proteins.

4. **Protein Microarrays:**
   - Immobilizes proteins on a surface and detects interactions using labeled probes.

5. **Fluorescence Resonance Energy Transfer (FRET):**
   - Measures energy transfer between two fluorophores attached to interacting proteins.

**Importance of PPIs:**
- **Cellular Function:** PPIs are involved in almost every cellular process, including signal transduction, immune responses, and cell cycle regulation.
- **Disease Mechanisms:** Abnormal PPIs are linked to diseases such as cancer, Alzheimer's, and viral infections.
- **Drug Targets:** Targeting specific PPIs can lead to the development of novel therapeutics.

**Databases and Resources:**
- **STRING:** Database of known and predicted PPIs.
- **BioGRID:** Repository of genetic and protein interactions.
- **IntAct:** Database for molecular interaction data.

**Visualization of PPIs:**
- **Network Diagrams:** Represent proteins as nodes and interactions as edges.
- **Cytoscape:** Software platform for visualizing complex networks.

**Challenges in PPI Research:**
- **False Positives/Negatives:** Experimental methods may produce inaccurate results.
- **Dynamic Nature:** PPIs can change under different cellular conditions.
- **Complexity:** High number of interactions makes analysis and interpretation difficult.
""",
        "tags": ["Proteomics", "Protein-Protein Interactions", "PPIs", "Yeast Two-Hybrid", "Co-IP"],
    },

    # Cell Biology
    {
        "id": 13,
        "category": "Biological Sciences",
        "title": "Cell Biology - Cell Structure and Organelles",
        "content": """
**Cell Structure and Organelles**

**Cell Types:**
1. **Prokaryotic Cells:**
   - Lack a nucleus and membrane-bound organelles.
   - Example: Bacteria and Archaea.

2. **Eukaryotic Cells:**
   - Possess a nucleus and various membrane-bound organelles.
   - Found in animals, plants, fungi, and protists.

**Key Organelles in Eukaryotic Cells:**
1. **Nucleus:**
   - Contains the cell's genetic material (DNA).
   - Enclosed by a nuclear envelope with nuclear pores.

2. **Mitochondria:**
   - Powerhouses of the cell, generating ATP through cellular respiration.
   - Double-membraned with their own DNA.

3. **Endoplasmic Reticulum (ER):**
   - **Rough ER:** Studded with ribosomes; involved in protein synthesis and folding.
   - **Smooth ER:** Lacks ribosomes; involved in lipid synthesis and detoxification.

4. **Golgi Apparatus:**
   - Modifies, sorts, and packages proteins and lipids for secretion or delivery to other organelles.

5. **Lysosomes:**
   - Contain digestive enzymes to break down waste materials and cellular debris.

6. **Peroxisomes:**
   - Break down fatty acids and amino acids; detoxify harmful substances.

7. **Ribosomes:**
   - Sites of protein synthesis.
   - Can be free-floating in the cytoplasm or bound to the rough ER.

8. **Cytoskeleton:**
   - Network of protein filaments (microtubules, actin filaments, intermediate filaments) providing structure, shape, and facilitating movement.

9. **Centrosomes and Centrioles:**
   - Involved in organizing microtubules during cell division.

10. **Vacuoles:**
    - Storage organelles; larger in plant cells (central vacuole) for maintaining turgor pressure.

11. **Chloroplasts (in Plant Cells):**
    - Sites of photosynthesis.
    - Contain chlorophyll and their own DNA.

12. **Cell Membrane:**
    - Phospholipid bilayer with embedded proteins.
    - Regulates the movement of substances in and out of the cell.

**Cell Wall (in Plant Cells):**
- Rigid outer layer providing structural support and protection.
- Composed mainly of cellulose.

**Cytoplasm:**
- Gel-like substance within the cell membrane containing organelles and cytosol.

**Nucleolus:**
- Located within the nucleus.
- Site of ribosomal RNA (rRNA) synthesis and ribosome assembly.

**Membrane Trafficking:**
- Process involving the movement of proteins and lipids between organelles via vesicles.
""",
        "tags": ["Cell Biology", "Organelles", "Eukaryotic Cells", "Prokaryotic Cells", "Cytoskeleton"],
    },
    {
        "id": 14,
        "category": "Biological Sciences",
        "title": "Cell Biology - Cell Signaling Pathways",
        "content": """
**Cell Signaling Pathways**

**Definition:**
- Complex systems of communication that govern basic cellular activities and coordinate cell actions.
- Enable cells to respond to their environment and maintain homeostasis.

**Key Components:**
1. **Signal Molecules (Ligands):**
   - Substances that initiate signaling cascades.
   - Examples: Hormones, growth factors, neurotransmitters.

2. **Receptors:**
   - Proteins on the cell surface or within the cell that bind to ligands.
   - **Types:**
     - **G-Protein Coupled Receptors (GPCRs):** Involved in many physiological processes.
     - **Receptor Tyrosine Kinases (RTKs):** Key in growth and differentiation.
     - **Intracellular Receptors:** Bind ligands that can cross the cell membrane (e.g., steroid hormones).

3. **Second Messengers:**
   - Small molecules that propagate the signal within the cell.
   - Examples: cAMP, Ca²⁺ ions, IP₃.

4. **Protein Kinases and Phosphatases:**
   - **Kinases:** Add phosphate groups to proteins, activating or deactivating them.
   - **Phosphatases:** Remove phosphate groups.

5. **Transcription Factors:**
   - Proteins that regulate gene expression in response to signaling.
   - Activate or repress the transcription of specific genes.

**Major Cell Signaling Pathways:**
1. **MAPK/ERK Pathway:**
   - Involved in cell growth, differentiation, and survival.
   - Triggered by growth factors binding to RTKs.

2. **PI3K/Akt Pathway:**
   - Regulates metabolism, growth, proliferation, and survival.
   - Activated by various growth factors and hormones.

3. **JAK/STAT Pathway:**
   - Mediates responses to cytokines and growth factors.
   - Involves Janus kinases (JAKs) and Signal Transducers and Activators of Transcription (STATs).

4. **Wnt/β-Catenin Pathway:**
   - Crucial for embryonic development and cell fate determination.
   - Dysregulation linked to cancers.

5. **Notch Signaling Pathway:**
   - Mediates cell-cell communication.
   - Important for cell differentiation processes.

6. **TGF-β Pathway:**
   - Regulates cell growth, differentiation, and immune responses.
   - Involved in development and wound healing.

**Examples of Cell Signaling in Action:**
- **Insulin Signaling:**
  - Insulin binds to its receptor, triggering the PI3K/Akt pathway.
  - Facilitates glucose uptake in cells.

- **Immune Response:**
  - Cytokines bind to receptors on immune cells, activating signaling pathways that regulate immune responses.

- **Neurotransmission:**
  - Neurotransmitters bind to receptors on neurons, initiating signaling cascades that influence neuronal activity.

**Dysregulation of Cell Signaling:**
- Can lead to diseases such as cancer, diabetes, and autoimmune disorders.
- Targeting specific components of signaling pathways is a strategy in drug development.
""",
        "tags": ["Cell Biology", "Cell Signaling", "Signaling Pathways", "GPCR", "RTK"],
    },
    {
        "id": 15,
        "category": "Biological Sciences",
        "title": "Cell Biology - Cell Cycle and Division",
        "content": """
**Cell Cycle and Division**

**Cell Cycle Overview:**
- **Definition:** Series of events that take place in a cell leading to its division and duplication.
- **Phases:**
  1. **Interphase:** Period of cell growth and DNA replication.
     - **G₁ Phase:** Cell growth and normal metabolic roles.
     - **S Phase:** DNA synthesis and replication.
     - **G₂ Phase:** Preparation for mitosis.
  2. **Mitosis (M Phase):** Division of the nucleus and its contents.
     - **Prophase:** Chromosomes condense; nuclear envelope breaks down.
     - **Metaphase:** Chromosomes align at the cell equator.
     - **Anaphase:** Sister chromatids are pulled apart to opposite poles.
     - **Telophase:** Nuclear envelopes reform around separated chromatids.
  3. **Cytokinesis:** Division of the cytoplasm, resulting in two daughter cells.

**Regulation of the Cell Cycle:**
- **Cyclins and Cyclin-Dependent Kinases (CDKs):**
  - Cyclins are proteins that regulate the activity of CDKs.
  - **Function:** Activate CDKs to trigger cell cycle transitions.
- **Checkpoints:**
  - **G₁ Checkpoint:** Ensures the cell is ready for DNA synthesis.
  - **G₂ Checkpoint:** Ensures DNA replication is complete and accurate.
  - **Metaphase Checkpoint:** Ensures all chromosomes are properly aligned before separation.

**Key Regulators:**
1. **p53 Protein:**
   - Acts as a tumor suppressor.
   - Initiates cell cycle arrest or apoptosis in response to DNA damage.

2. **Retinoblastoma Protein (Rb):**
   - Regulates the G₁/S transition.
   - Phosphorylated by Cyclin D/CDK4/6, allowing progression to S phase.

3. **Anaphase-Promoting Complex/Cyclosome (APC/C):**
   - E3 ubiquitin ligase that targets proteins for degradation, allowing mitotic exit.

**Types of Cell Division:**
1. **Mitosis:**
   - Results in two genetically identical diploid daughter cells.
   - Responsible for growth, repair, and asexual reproduction.

2. **Meiosis:**
   - Specialized form of cell division producing four genetically diverse haploid gametes.
   - Involves two successive divisions: Meiosis I and Meiosis II.
   - Essential for sexual reproduction and genetic diversity.

**Cytokinesis Mechanisms:**
- **Animal Cells:**
  - **Cleavage Furrow:** Actin and myosin filaments contract, pinching the cell into two.
- **Plant Cells:**
  - **Cell Plate Formation:** Vesicles coalesce to form a new cell wall separating the daughter cells.

**Cell Cycle Dysregulation:**
- **Cancer:** Uncontrolled cell division due to mutations in cell cycle regulators.
- **Genetic Disorders:** Errors in cell division can lead to conditions like Down syndrome.

**Techniques to Study the Cell Cycle:**
- **Flow Cytometry:** Analyzes cell cycle distribution by measuring DNA content.
- **Live-Cell Imaging:** Observes cell division in real-time.
- **Fluorescence Microscopy:** Visualizes specific cell cycle proteins and structures.
""",
        "tags": ["Cell Biology", "Cell Cycle", "Mitosis", "Cytokinesis", "Cell Division"],
    },

    # Evolution
    {
        "id": 16,
        "category": "Biological Sciences",
        "title": "Evolution - Natural Selection and Genetic Drift",
        "content": """
**Natural Selection and Genetic Drift**

**Natural Selection:**
- **Definition:** Process where organisms better adapted to their environment tend to survive and produce more offspring.
- **Types:**
  1. **Directional Selection:** Favors one extreme phenotype.
  2. **Stabilizing Selection:** Favors intermediate phenotypes.
  3. **Disruptive Selection:** Favors both extreme phenotypes.

- **Mechanism:**
  - Variation in traits exists within a population.
  - Environmental pressures select for advantageous traits.
  - Over generations, advantageous traits become more common.

- **Example:**
  - Peppered Moth during the Industrial Revolution: Dark-colored moths became more common in polluted areas due to better camouflage against predators.

**Genetic Drift:**
- **Definition:** Random changes in allele frequencies in a population, especially pronounced in small populations.
- **Causes:**
  - **Bottleneck Effect:** Reduction in population size due to events like natural disasters.
  - **Founder Effect:** New population started by a small number of individuals.

- **Effects:**
  - Can lead to the loss of genetic variation.
  - Alleles may become fixed or lost purely by chance.

- **Example:**
  - In a small population of birds, a random event like a storm may eliminate most individuals, drastically changing allele frequencies.

**Comparison Between Natural Selection and Genetic Drift:**
| Aspect               | Natural Selection                           | Genetic Drift                        |
|----------------------|---------------------------------------------|--------------------------------------|
| **Cause**            | Differential survival and reproduction      | Random sampling of alleles           |
| **Effect on Allele Frequency** | Directional based on trait advantage | Random, unpredictable changes        |
| **Population Size** | Effective in large populations              | More pronounced in small populations |
| **Role in Evolution** | Leads to adaptation and increased fitness | Can lead to loss of genetic variation |

**Other Evolutionary Forces:**
- **Mutation:** Introduction of new genetic variations.
- **Gene Flow:** Movement of alleles between populations.

**Importance in Evolution:**
- **Natural Selection:** Drives adaptive evolution, leading to traits that enhance survival and reproduction.
- **Genetic Drift:** Contributes to genetic diversity and can influence evolutionary paths, especially in isolated or small populations.

**Real-World Implications:**
- **Conservation Biology:** Understanding genetic drift is crucial for managing endangered species.
- **Agriculture:** Natural selection principles are applied in selective breeding programs.
- **Medicine:** Insights into evolutionary forces aid in understanding antibiotic resistance.

**Conclusion:**
Both natural selection and genetic drift are fundamental mechanisms of evolution, shaping the genetic makeup of populations over time through different processes.
""",
        "tags": ["Evolution", "Natural Selection", "Genetic Drift", "Adaptive Evolution", "Population Genetics"],
    },
    {
        "id": 17,
        "category": "Biological Sciences",
        "title": "Evolution - Molecular Evolution",
        "content": """
**Molecular Evolution**

**Definition:**
- Study of evolutionary changes at the molecular level, focusing on DNA, RNA, and protein sequences.
- Examines the mechanisms driving genetic variation and divergence among species.

**Key Concepts:**
1. **Mutation:**
   - Changes in the nucleotide sequence.
   - Sources of genetic variation.

2. **Selection Pressures:**
   - **Positive Selection:** Favours advantageous mutations.
   - **Purifying Selection:** Removes deleterious mutations.
   - **Neutral Selection:** No effect on fitness; molecular drift dominates.

3. **Genetic Code and Codon Usage:**
   - Evolutionary conservation and variability in codon preferences.
   - Codon bias can affect gene expression levels.

4. **Molecular Clocks:**
   - Estimate the time of evolutionary divergence based on genetic mutations.
   - Assumes a constant rate of mutation over time.

5. **Homologous Genes:**
   - **Orthologs:** Genes in different species evolved from a common ancestral gene.
   - **Paralogs:** Genes related by duplication within a genome.

6. **Gene Family Evolution:**
   - Expansion and contraction of gene families through duplication and loss.
   - Contributes to functional diversification.

**Mechanisms of Molecular Evolution:**
1. **Point Mutations:**
   - Single nucleotide changes affecting protein function or regulation.
   
2. **Insertions and Deletions (Indels):**
   - Addition or loss of nucleotides altering the reading frame or gene structure.
   
3. **Recombination:**
   - Exchange of genetic material between DNA molecules.
   - Increases genetic diversity.

4. **Horizontal Gene Transfer:**
   - Transfer of genes between unrelated species.
   - Common in prokaryotes, contributing to antibiotic resistance.

**Molecular Phylogenetics:**
- Reconstruction of evolutionary relationships using molecular data.
- Utilizes DNA, RNA, or protein sequences to build phylogenetic trees.

**Applications:**
- **Understanding Evolutionary Relationships:** Clarifying the evolutionary history of species.
- **Identifying Functional Regions:** Conserved sequences often indicate essential functions.
- **Disease Research:** Tracing the evolution of pathogens and resistance genes.
- **Biotechnology:** Designing enzymes and proteins with desired properties based on evolutionary insights.

**Tools and Databases:**
- **BLAST:** Aligns sequences to identify homologs.
- **PhyML, RAxML:** Software for building phylogenetic trees.
- **PAML:** Software for phylogenetic analysis by maximum likelihood.

**Challenges:**
- **Rate Variability:** Mutation rates can vary among lineages, complicating molecular clock estimates.
- **Sequence Alignment:** Accurate alignment is crucial for reliable phylogenetic inference.
- **Horizontal Gene Transfer:** Can obscure true evolutionary relationships.

**Conclusion:**
Molecular evolution provides deep insights into the genetic basis of diversity, adaptation, and the intricate history of life on Earth by analyzing the changes and dynamics at the molecular level.
""",
        "tags": ["Molecular Evolution", "Phylogenetics", "Genetic Variation", "Evolutionary Biology"],
    },
    {
        "id": 18,
        "category": "Biological Sciences",
        "title": "Evolution - Phylogenetics",
        "content": """
**Phylogenetics**

**Definition:**
- Branch of biology that deals with the evolutionary development and diversification of a species or group of organisms.
- Constructs phylogenetic trees (cladograms) to represent evolutionary relationships.

**Key Concepts:**
1. **Phylogenetic Trees:**
   - **Structure:** Nodes represent common ancestors, and branches represent evolutionary lineages.
   - **Types:**
     - **Cladograms:** Show relationships based on shared derived characteristics.
     - **Phylograms:** Branch lengths represent the amount of evolutionary change.
     - **Ultrametrical Trees:** Branch lengths represent time.

2. **Cladistics:**
   - Method of classification based on common ancestry.
   - Uses shared derived traits (synapomorphies) to infer relationships.

3. **Parsimony:**
   - Principle that the simplest explanation (fewest evolutionary changes) is preferred.
   - Used in tree construction algorithms like Maximum Parsimony.

4. **Maximum Likelihood:**
   - Statistical method that finds the tree most likely to have produced the observed data.
   - Accounts for different rates of evolution.

5. **Bayesian Inference:**
   - Probabilistic approach incorporating prior knowledge.
   - Provides a posterior probability distribution of trees.

6. **Molecular Phylogenetics:**
   - Utilizes molecular data (DNA, RNA, protein sequences) to reconstruct evolutionary relationships.
   - Often more precise than morphological data.

**Tools and Software:**
- **MEGA (Molecular Evolutionary Genetics Analysis):** Comprehensive tool for phylogenetic analysis.
- **RAxML (Randomized Axelerated Maximum Likelihood):** Software for large-scale phylogenetic analyses.
- **MrBayes:** Bayesian inference tool for phylogenetics.
- **PhyML:** Software for maximum likelihood phylogenetic tree estimation.

**Applications:**
- **Understanding Evolutionary
