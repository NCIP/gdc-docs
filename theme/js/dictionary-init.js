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

        if (currentView.indexOf('entity') < 0) {
          dictionaryPreamble.hide();

          if (previousView !== currentView || window.location.hash.indexOf('_top') >= 0) {
            body.scrollTop(0);
          }
        }
        else {
          dictionaryPreamble.show();
        }

        previousView = currentView;

      }
    };

    window.$gdcApp.dictionaryViewer = new Dictionary(dictionaryContainer, dictionaryOptions);

  }

  _init();
};

window.onunload = function() {
  window.$gdcApp.dictionaryViewer.destroy();
};