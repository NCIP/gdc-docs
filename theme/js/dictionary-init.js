/*
window.onload = function() {
  'use strict';

  window.gdcApp = window.gdcApp || {};

  function _init() {
    var dictionaryContainer = document.getElementById('dictionary-app-container');

    if (! dictionaryContainer ) {
      return;
    }

    var dictionaryOptions = {
      dataSourceBaseHost: 'http://localhost:8080'
    };

    window.$gdcApp.dictionaryViewer = new Dictionary(dictionaryContainer, dictionaryOptions);

  }

  _init();
};

window.onunload = function() {
  window.$gdcApp.dictionaryViewer.destroy();
};
*/