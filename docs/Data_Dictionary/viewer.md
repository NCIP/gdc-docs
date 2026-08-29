---
hide:
    - navigation
    - toc
---

# Data Dictionary Viewer

<link rel="stylesheet" href="dictionary.css">
<script src="https://cdn.jsdelivr.net/npm/lodash@4.17.21/lodash.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js" integrity="sha512-vc58qvvBdrDR4etbxMdlTt4GBQk1qjvyORR2nrsPsFPyrs+/u5c3+1Ct6upOgdZoIl7eq6k3a1UPDSNAQi/32A==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>

<script src="dictionary-init.js"></script>
<script src="dictionary.js"></script>
<script src="dictionary-views.js"></script>
<p id="dictionary-preamble" style="margin: 1rem auto 2rem auto; width: 95%;">
The GDC data dictionary viewer is a user-friendly interface for accessing the <a href="../">GDC Data Dictionary</a>.
</p>
<div id="dictionary-loading-icon" class="loadingContainer">
    <div class="spinParticleContainer">
        <div class="particle red"></div>
        <div class="particle grey other-particle"></div>
        <div class="particle blue other-other-particle"></div>
    </div>
  <div>
  Loading Dictionary Data...
  </div>
</div>
<div id="dictionary-app-container"></div>
