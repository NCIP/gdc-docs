
<!-- GDCMVS APP -->

<h1 class="hide">GDC Metadata Validation Services</h1>
<div class="search-box">
    <div class="search-bar">
      <div class="input-group search-bar__input">
        <input id="keywords" type="text" class="form-control" aria-label="keywords" placeholder="Please enter keyword">
        <a href="#" id="searchclear" class="search-bar__clear" aria-label="clear search bar" style="display: none;"><i class="fa fa-times"></i></a>
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
        <a href="https://cdebrowser.nci.nih.gov/" class="ref-box__link" target="_blank">Search in caDSR</a>
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

<div class="version-content">GDC Dictionary Version - Jibboo</div>
<div id="alert-error" class="alert alert__error alert-info alert-dismissible" role="alert" style="display: none;">Error: Undefined</div>

<!-- END GDCMVS APP -->
