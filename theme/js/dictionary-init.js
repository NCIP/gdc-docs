window.onload = function() {
  'use strict';

  window.gdcApp = window.gdcApp || {};

  // For all multi-value options identified by the dictionary
  // create a simple accordion component that will allow us to
  // hide/show more values
  function _initMultiValueAccordions(numberOfShowingItems) {

    var accordions = jQuery('.bullets');

    if (accordions.length === 0) {
      return;
    }

    accordions.parent().append('<div class="accordion-more-less-toggle" />');

    var bullets = accordions.find('li li');

    if (bullets.length === 0) {
      return;
    }

    var liHeight = bullets.outerHeight(),
      valueHeight = numberOfShowingItems * liHeight;

    accordions.each(function () {
      var accordion = jQuery(this);

      if (accordion.outerHeight() > valueHeight) {
        accordion.css({overflow: 'hidden', height: valueHeight});
        accordion.parent().find('.accordion-more-less-toggle')
          .css({marginLeft: '2.2rem'})
          .append('<a href="javascript:void(0)" class="more-less-vals" title="Click to see more values."><i class="fa fa-angle-down"></i> More Values</a>');
      }
    });

    jQuery('.more-less-vals').click(function () {
      var moreLessToggle = jQuery(this),
        bullets = moreLessToggle.parent().parent().find('.bullets'),
        title = 'Click to see less values.',
        html = '<i class="fa fa-angle-up"></i> Less Values</a>',
        css = {overflow: 'visible', height: ''};

      if (moreLessToggle.hasClass('more')) {
        css = {overflow: 'hidden', height: valueHeight};
        html = '<i class="fa fa-angle-down"></i> More Values</a>';
      }

      bullets.css(css);
      moreLessToggle.attr('title', title).html(html);

      moreLessToggle.toggleClass('more');
    });

  }

  //////////////////////////////////////////////////////
  // Init function for our dictionary
  //////////////////////////////////////////////////////
  function _init() {
    var dictionaryContainer = document.getElementById('dictionary-app-container');

    if (! dictionaryContainer ) {
      return;
    }

    var dictionaryPreamble = jQuery('#dictionary-preamble'),
        body = jQuery('body'),
        previousView = null;

    var dictionaryOptions = {
      //dataSourceBaseHost: 'http://localhost:8080',
      afterRenderFn: function(dictionary) {

        var currentView = dictionary.getCurrentViewName();
        var shouldScrollToTop = window.location.hash.indexOf('_top') >= 0;

        if (currentView !== dictionary.CONSTANTS.VIEWS.TABLE.ENTITY_LIST) {
          dictionaryPreamble.hide();

          if (previousView !== currentView) {
            body.scrollTop(0);
          }
        }
        else {
          dictionaryPreamble.show();
        }

        if (shouldScrollToTop) {
          body.scrollTop(0);
        }

        previousView = currentView;

        _initMultiValueAccordions(6); // 5 values + 1 header

      }
    };

    window.$gdcApp.dictionaryViewer = new Dictionary(dictionaryContainer, dictionaryOptions);

  }

  _init();
};

window.onunload = function() {
  window.$gdcApp.dictionaryViewer.destroy();
};