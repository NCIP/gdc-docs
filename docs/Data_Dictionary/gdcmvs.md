
<!-- GDCMVS APP -->

<h1 class="hide">GDC Metadata Validation Services</h1>
<div class="search-box">
    <div class="search-bar">
      <span id="suggestWidth" class="suggest__width"></span>
      <div class="input-group search-bar__input">
        <input id="keywords" type="text" class="form-control" aria-label="keywords" placeholder="Please enter keyword" autocomplete="off">
        <div id="search-bar-options" class="search-bar__options dropdown" style="display: none;">
          <a href="#" data-toggle="dropdown" class="dropdown-toggle search-bar__option" aria-label="boolean operators"><i class="fa fa-ellipsis-h"></i></a>
          <ul class="dropdown-menu search-bar__dropdown">
              <li><a class="search-bar__boolean" data-boolean="AND" href="#">AND</a></li>
              <li><a class="search-bar__boolean" data-boolean="OR" href="#">OR</a></li>
              <li><a class="search-bar__boolean" data-boolean="NOT" href="#">NOT</a></li>
          </ul>
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
        <a href="https://cdebrowser.nci.nih.gov/cdebrowserClient/cdeBrowser.html#/search?programArea=0&contextId=2C8BAF10-7E19-B797-E050-BB89AD43619C"
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
