# GDC Encyclopedia
<link rel="stylesheet" type="text/css" href="encyclopedia.css">
<div style="display: flex; flex-direction: column"; class="encyclopedia-container">
    <div style="display: flex; flex direction column;" class = "encyclopedia-main">
        <div style="display: flex; flex-direction: column;" class="encyclopedia-menu">
            <div class="side-box">
                <span class="side-box-title">Search the Encyclopedia<br>
                <input type="text" style="margin-bottom: 12px;" class="side-box-input encyclopedia-search">
                </span>
                <span class="side-box-title">Browse the terms</span>
                <br>
                <div style="line-height: 34px; margin-top: 10px; margin-bottom: 7px">
                    <button class='term-button term-link' id="A-button">A</button>
                    <button class='term-button term-link' id="B-button">B</button>
                    <button class='term-button term-link' id="C-button">C</button>
                    <button class='term-button term-link' id="D-button">D</button>
                    <button class='term-button term-link' id="E-button">E</button>
                    <button class='term-button term-link' id="F-button">F</button>
                    <button class='term-button term-link' id="G-button">G</button>
                    <button class='term-button term-link' id="H-button">H</button>
                    <button class='term-button term-link' id="I-button">I</button>
                    <button class='term-button term-link' id="J-button">J</button>
                    <button class='term-button term-link' id="K-button">K</button>
                    <button class='term-button term-link' id="L-button">L</button>
                    <button class='term-button term-link' id="M-button">M</button>
                    <button class='term-button term-link' id="N-button">N</button>
                    <button class='term-button term-link' id="O-button">O</button>
                    <button class='term-button term-link' id="P-button">P</button>
                    <button class='term-button term-link' id="Q-button">Q</button>
                    <button class='term-button term-link' id="R-button">R</button>
                    <button class='term-button term-link' id="S-button">S</button>
                    <button class='term-button term-link' id="T-button">T</button>
                    <button class='term-button term-link' id="U-button">U</button>
                    <button class='term-button term-link' id="V-button">V</button>
                    <button class='term-button term-link' id="W-button">W</button>
                    <button class='term-button term-link' id="X-button">X</button>
                    <button class='term-button term-link' id="Y-button">Y</button>
                    <button class='term-button term-link' id="Z-button">Z</button>
                    <button class='term-button term-link' id="NUM-button">#</button>
                </div>
                <a style="font-weight: bold" class="side-box-text" id="view-all-link" href="#">View All <i class="fa fa-chevron-right"></i></a> 
            </div>
            <div class="side-box">
                <span class="side-box-title" style="padding-right: 108px; border-bottom: 1px solid #e5e5e5;">Suggest a Topic</span>
                <br>
                <br>
                <span class="side-box-text" style="line-height: 10px;">Our goal is to maintain this resource as an important tool for understanding cancer biology. If there is
                    a topic you are unsure about, please
                    <a href="" style="color: #6385a2"> contact us</a> or suggest the topic below.
                </span> 
                <br>
                <br>
                <span class="side-box-text" style="font-weight: bold">Suggest a topic:</span>
                <input type="text" class="side-box-input topic-suggest">
                <button class="topic-submit-button">Submit</button>
            </div>
        </div>
        <div class="encyclopedia-content" style="margin-left: 10px;">
            <div>
                The GDC Encyclopedia is an informational tool that makes GDC applications easier to use. It helps researchers, data submitters, developers and clinicians quickly find specific information without needing to browse for it.
            </div>    
            <div id = "encyc-table" style="color: white">
                <header id = "encyc-table-header">
                    <h3>Common Topics for Researchers</h3>
                </header>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Controlled access</a>
                    </span>
                    <span>
                        <a href="">dbGaP</a>
                    </span>
                    <span>
                        <a href="">Open access</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Data Browser</a>
                    </span>
                    <span>
                        <a href="">Experiment</a>
                    </span>
                    <span>
                        <a href="">Tissue Source Site (TSS)</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Data Matrix</a>
                    </span>
                    <span>
                        <a href="">Investigation Description Format (IDF)</a>
                    </span>
                    <span>
                        <a href="">Underexpressed</a>
                    </span>
                </div>
            </div>
            <div id = "encyc-table" style="color: white">
                <header id = "encyc-table-header">
                    <h3>Common Topics for Data Submitters</h3>
                </header>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Data Submission Center</a>
                    </span>
                    <span>
                        <a href="">File Type</a>
                    </span>
                    <span>
                        <a href="">Pipeline</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">dbGaP</a>
                    </span>
                    <span>
                        <a href="">Latest Archive Report</a>
                    </span>
                    <span>
                        <a href="">Quality management system</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">File Format</a>
                    </span>
                    <span>
                        <a href="">Normalization</a>
                    </span>
                    <span>
                        <a href="">TCGA barcode</a>
                    </span>
                </div>
            </div>
            <div id = "encyc-table" style="color: white">
                <header id = "encyc-table-header">
                    <h3>Common Topics for Data Developers</h3>
                </header>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Facets</a>
                    </span>
                    <span>
                        <a href="">Gzip</a>
                    </span>
                    <span>
                        <a href="">Normalization</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">FASTA Format</a>
                    </span>
                    <span>
                        <a href="">HTTP data access</a>
                    </span>
                    <span>
                        <a href="">QCLive</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Genome Data Analysis Center (GDAC)</a>
                    </span>
                    <span>
                        <a href="">MD Anderson GDAC MBatch</a>
                    </span>
                    <span>
                        <a href="">Universally Unique Identifier</a>
                    </span>
                </div>
            </div>
            <div id = "encyc-table" style="color: white">
                <header id = "encyc-table-header">
                    <h3>Common Topics for Clinicians</h3>
                </header>
                <div id="encyc-table-content">
                    <span>
                        <a href="">DNA methylation</a>
                    </span>
                    <span>
                        <a href="">HUGO gene symbol</a>
                    </span>
                    <span>
                        <a href="">Metadata</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Genome Characterization Center</a>
                    </span>
                    <span>
                        <a href="">Loss of heterozygosity (LOH)</a>
                    </span>
                    <span>
                        <a href="">Protected Health Information</a>
                    </span>
                </div>
                <div id="encyc-table-content">
                    <span>
                        <a href="">Gene expression data (array based)</a>
                    </span>
                    <span>
                        <a href="">Mutation Annotation Format (MAF)</a>
                    </span>
                    <span>
                        <a href="">Tumor mutation samples</a>
                    </span>
                </div>
            </div>
        </div>
    </div>
</div>