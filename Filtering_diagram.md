
flowchart TD
    RAW["00 — Raw FreeBayes variants<br/>24,037,651 sites | 263 individuals"]
    BC["01 — Initial BCFtools<br/>17,907,986 SNPs"]
    VC["02 — Initial VCFtools<br/>8,933,146 SNPs"]
    IND["03 — Individual missingness<br/>All 263 assessed"]
    REM["04a — Remove 5 individuals<br/>258 retained | 8,933,146 SNPs"]
    SITE["04b — Site completeness ≥70%<br/>1,334,122 SNPs"]
    MAC["05 — Reapply MAC ≥2<br/>1,332,396 SNPs"]
    DEP["06 — Mean site depth assessment<br/>No SNPs removed"]
    SPLIT{"07 — Global MAF<br/>Two datasets"}

    RAW --> BC --> VC --> IND --> REM --> SITE --> MAC --> DEP --> SPLIT

    IND -.-> IOLD["Earlier approach<br/>Keep all 263 individuals<br/>for initial diagnostics"]
    IOLD -.-> M263["Original site completeness assessment<br/>263 individuals"]
    M263 -.-> THR["Alternative site-presence thresholds<br/>95%: 41,854 SNPs<br/>90%: 119,552<br/>80%: 613,206<br/>70%: 1,256,391<br/>50%: 2,388,252"]
    THR -.-> REASSESS["Reassessment<br/>Remove 5 high-missingness individuals<br/>Recalculate site completeness"]
    REASSESS -.-> REM

    DEP -.-> DNOTE["Decision<br/>No additional minimum mean-depth filter<br/>Maximum mean depth evaluated later"]

    SPLIT --> N7["Neutral — MAF ≥0.05<br/>536,457 SNPs"]
    SPLIT --> A7["Adaptive — MAF ≥0.01<br/>879,199 SNPs"]

    SPLIT -.-> MALT["Original MAF assessment<br/>0.01: 879,199<br/>0.02: 687,360<br/>0.03: 624,857<br/>0.05: 536,457"]
    MALT -.-> MDEC["Original decision: MAF ≥0.01 only<br/>Later revised to two parallel datasets"]
    MDEC -.-> SPLIT

    N7 --> N8["08 — Mean depth ≤500×<br/>536,448 SNPs"]
    A7 --> A8["08 — Mean depth ≤500×<br/>879,186 SNPs"]

    N8 --> N9["09 — Single-base alleles<br/>501,299 SNPs"]
    A8 --> A9["09 — Single-base alleles<br/>839,809 SNPs"]

    N9 --> N10["10 — Revised AB filtering<br/>0.30 &lt; AB &lt; 0.75 or AB &lt; 0.01<br/>Count pending"]
    A9 --> A10["10 — Revised AB filtering<br/>0.30 &lt; AB &lt; 0.75 or AB &lt; 0.01<br/>Count pending"]

    A9 -.-> ABOLD["Original AB filter<br/>0.25 &lt; AB &lt; 0.75 or AB &lt; 0.01<br/>628,402 SNPs"]
    ABOLD -.-> MQOLD["Original mapping-quality filter<br/>595,806 SNPs"]
    ABOLD -.-> ABREV["Advisor review<br/>Clear break at AB = 0.30<br/>Lower threshold revised"]
    ABREV -.-> A10
    ABREV -.-> N10

    N10 --> N11["11 — Mapping quality<br/>Count pending"]
    A10 --> A11["11 — Mapping quality<br/>Count pending"]

    N11 --> N12["12 — Heterozygosity assessment<br/>Pending"]
    A11 --> A12["12 — Heterozygosity assessment<br/>Pending"]

    classDef shared fill:#e5effc,stroke:#4777b2,color:#173b67
    classDef neutral fill:#e4f3e8,stroke:#3d8655,color:#225034
    classDef adaptive fill:#e9e6fa,stroke:#7863b2,color:#3c306b
    classDef history fill:#fff0d8,stroke:#bd873e,color:#6d481a
    classDef pending fill:#f2f2f2,stroke:#888,color:#444

    class RAW,BC,VC,IND,REM,SITE,MAC,DEP shared
    class N7,N8,N9 neutral
    class A7,A8,A9 adaptive
    class IOLD,M263,THR,REASSESS,DNOTE,MALT,MDEC,ABOLD,MQOLD,ABREV history
    class N10,A10,N11,A11,N12,A12 pending
