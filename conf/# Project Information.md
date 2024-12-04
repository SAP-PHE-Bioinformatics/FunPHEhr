# Project Information

**Project Name:**  
**Principal Investigator:**  
**Contact Information:**  
**Funding Source:**  
**Project Duration:**  
**Data Collection Date:**  

# HPC System Information

**HPC Facility Name:**  
**System Configuration:**  
  - Number of Nodes: 3
  - Nodes: 4
  - slurm workload manager v23.11.3
  - SAN storage HDD: 96TB (vg_me4012-lv_phe)
  - **Nodes:**
  
   -frgeneseq03: 
    - Control node as well
    - CPUs 96
    - RAM: 386 GB
    - Storage Details:  
      - local /scratch ssd : 11TB
    
    -frgeneseq-control : 
    - Control node as well
    - CPUs 64
    - RAM: 192 GB
    - Storage Details:  
      - local /scratch ssd : 200GB

    -frgeneseqgpu (Goliath) : 
    - Control node as well
    - CPUs 64
    - RAM: 128  GB
    - GRES GPU: 2x Nvidia 
    - Storage Details:  
      - local /scratch ssd : 3.6TB
      - local secondary /data ssd raid : 7.4TB
      - local cold storage /slowNsteady HDD 15TB

    -David : 
    - Control node as well
    - CPUs 48
    - RAM: 128  GB
    - GRES GPU: 1x Nvidia 3090Ti
    - Storage Details:  
      - local /scratch ssd : 3.6TB
      - local secondary /data ssd raid : 7.4TB
      - local cold storage /slowNsteady HDD 15TB

  
**Access Method:** (e.g., SSH, Remote Desktop)

# Pipeline and Workflow Details

**Pipeline Name:** Phesiqcal
**Version:**  
**Purpose of the Pipeline:** Bacterial de-Novo assembly, QC assessment and Annotation
**Major Software and Tools Used:**  
  - Tool Name:  
  - Version:  
  - Configuration:  
**Custom Scripts/Modifications:**  

**Pipeline Name:** PHET
**Version:**  
**Purpose of the Pipeline:** Bacterial Typing workflow, covers large spectrum of bacteria and runs speicfic typing tools for each genus. 
**Major Software and Tools Used:**  
  - Tool Name:  
  - Version:  
  - Configuration:  
**Custom Scripts/Modifications:**  

**Pipeline Name:** Phesymphony
**Version:**  
**Purpose of the Pipeline:** Bacterial SNp and analysis and Phylogenetics clustering
**Major Software and Tools Used:**  
  - Tool Name:  
  - Version:  
  - Configuration:  
**Custom Scripts/Modifications:** 



# Workflow Execution Metrics

**Average Number of Samples per Run:**  
**Frequency of Runs:**  
  - Times Run per Week:  
  - Average Duration of Each Run:  

# Resource Usage

**CPU Hours Utilized:**  
**Memory Utilized (GB):**  
**Storage Utilized (TB):**  
**GPU Hours Utilized (if applicable):**  

# Input Data

**Data Source:**  
**Data Type:**  
**Data Volume:**  
**Data Format:**  

# Output Data

**Primary Results:**  
**Data Volume:**  
**Data Format:**  
**Storage Location:**  

# Analysis and Interpretation

**Summary of Findings:**  
**Data Analysis Tools Used:**  
  - Tool Name:  
  - Version:  
  - Purpose:  
**Challenges Faced:**  
**Adjustments to Pipeline:**  

# Compliance and Security

**Data Privacy Compliance (e.g., HIPAA, GDPR):**  
**Security Measures:**  
**Data Retention Policy:**  

# Notes and Additional Comments
