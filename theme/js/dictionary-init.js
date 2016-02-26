window.onload = function() {
  'use strict';

  window.gdcApp = window.gdcApp || {};

  function _init() {
    var dictionaryContainer = document.getElementById('dictionary-app-container');

    if (! dictionaryContainer ) {
      return;
    }

    var dictionaryPreamble = $('#dictionary-preamble'),
        body = $('body'),
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

        var accordions = jQuery('.bullets');

        if (accordions.length) {
          accordions.parent().append('<div class="accordion-more-less-toggle" />');

          var bullets = accordions.find('li li');

          if (bullets.length === 0) {
            return;
          }

          var liHeight = bullets.outerHeight(),
              valueHeight = 6 * liHeight;

          accordions.each(function() {
            var $this = $(this);

            if($this.outerHeight() > valueHeight) {
              $this.css({overflow: 'hidden', height: valueHeight});
              $this.parent().find('.accordion-more-less-toggle')
                .css({marginLeft: '2.2rem'})
                .append('<a href="javascript:void(0)" class="more-less-vals" title="Click to see more values."><i class="fa fa-angle-down"></i> More Values</a>');
            }
          })

          jQuery('.more-less-vals').click(function(){
            var $this = $(this);



            if($this.hasClass('more')){
              $this.parent().parent().find('.bullets').css({overflow: 'hidden', height: valueHeight});
              $this.attr('title', 'Click to see more values.').html('<i class="fa fa-angle-down"></i> More Values</a>');
            }
            else {
              $this.parent().parent().find('.bullets').css({overflow: 'visible', height: ''});
              $this.attr('title', 'Click to see less values.').html('<i class="fa fa-angle-up"></i> Less Values</a>');
            }

            $this.toggleClass('more');
          });
        }



      }
    };

    window.$gdcApp.dictionaryViewer = new Dictionary(dictionaryContainer, dictionaryOptions);

  }

  _init();
};

window.onunload = function() {
  window.$gdcApp.dictionaryViewer.destroy();
};