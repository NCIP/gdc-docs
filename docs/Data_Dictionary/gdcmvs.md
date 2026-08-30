---
hide:
    - navigation
    - toc
---

<script src="https://code.jquery.com/jquery-3.7.1.min.js" integrity="sha256-/JqT3SQfawRcv/BIHPThkBvs0OEvtFFmqPF/lYI/Cxo=" crossorigin="anonymous"></script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.bundle.min.js" integrity="sha384-YvpcrYf0tY3lHB60NNkmXc5s9fDVZLESaAA55NDzOxhy9GkcIdslK1eN7N6jIeHz" crossorigin="anonymous"></script>
<script type="text/javascript" src="lib/js/jquery-ui.min.js"></script>
<script type="text/javascript" src="lib/js/pagination.min.js"></script>
<script type="text/javascript" src="lib/js/jsrender.min.js"></script>
<!-- boostrap used by bundle.js for tooltips -->
<script type="text/javascript" src="dist/bundle.js"></script>
<link rel="stylesheet" href="lib/css/jquery-ui.min.css">
<link rel="stylesheet" href="lib/css/pagination.min.css">
<link rel="stylesheet" href="dist/styles.css">


<!-- GDCMVS APP -->

<h1 class="hide">GDC Metadata Validation Services</h1>
<div class="search-box">
    <div class="search-bar">
      <span id="suggestWidth" class="suggest__width"></span>
      <div class="input-group search-bar__group">
        <input id="keywords" type="text" class="form-control search-bar__input" aria-label="keywords" placeholder="Please enter keyword" autocomplete="off">
        <div id="search-bar-options" class="search-bar__options dropdown" style="display: none;">
          <a href="#" data-toggle="dropdown" class="dropdown-toggle search-bar__option" aria-label="boolean operators"><i class="fa fa-ellipsis-h"></i></a>
          <div class="dropdown-menu search-bar__dropdown">
              <a class="search-bar__boolean" data-boolean="AND" href="#">AND</a>
              <a class="search-bar__boolean" data-boolean="OR" href="#">OR</a>
              <a class="search-bar__boolean" data-boolean="NOT" href="#">NOT</a>
          </div>
          <a href="#" id="searchclear" class="search-bar__option" aria-label="clear search bar"><i class="fa fa-times"></i></a>
        </div>
        <span class="input-group-btn">
          <button id="search" class="btn search-bar__btn" type="button">Search</button>
        </span>
      </div>
      <div class="suggest">
        <div id="suggestBox" class="suggest__listbox"></div>
      </div>
    </div>
    <div class="search-options">
      <div class="checkbox">
        <label class="checkbox__label checkbox__label--padding">
          <input id="i_ematch" class="checkbox__input" type="checkbox" value="" tabindex="0">
          <span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Exact match
        </label>
        <label class="checkbox__label">
          <input id="i_desc" class="checkbox__input" type="checkbox" value="" tabindex="0">
          <span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Property description
        </label>
        <label class="checkbox__label">
          <input id="i_syn" class="checkbox__input" type="checkbox" value="" tabindex="0">
          <span class="checkbox__btn"><i class="checkbox__icon fa fa-check"></i></span> Synonyms
        </label>
      </div>
      <div class="ref-box">
        <a href="https://ncit.nci.nih.gov/" class="ref-box__link" target="_blank">Search in NCIt</a>
        <a href="https://cadsr.cancer.gov/onedata/Home.jsp"
          class="ref-box__link" target="_blank">Search in caDSR</a>
      </div>
    </div>
</div>

<div id="gdc-loading-icon" class="loadingContainer" style="display: none;">
  <div class="spinParticleContainer">
      <div class="particle red"></div>
      <div class="particle grey other-particle"></div>
      <div class="particle blue other-other-particle"></div>
  </div>
  <div>Loading GDC Data...</div>
</div>

<div id="root"></div>

<div id="info-content" class="info-content">
    <div id="unofficial-term"></div>
    <div id="version-content" class="version-content">GDC Dictionary Version</div>
</div>
<div id="alert-error" class="alert alert__error alert-info" role="alert">Error: Undefined</div>
<div id="alert-suggest" class="alert alert__suggest alert-info" role="alert">
  Coming soon: users will be able to suggest new vocabulary terms to the GDC via this link.
</div>

<!-- END GDCMVS APP -->
