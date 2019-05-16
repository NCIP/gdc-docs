$(function () {
  var $inputBox = $('#hp-search__input');
  var searchInput = document.getElementById('hp-search__input');
  var $submit = $('#hp-search__submit');
  var $resultsContainer = $('#hp-search__results-container');
  var $results = $('#hp-search__results');
  var $searchContentBody = $('#hp-search__body');
  var _VALID_QUERY_LENGTH = 3;
  var _isSearchActive = false;

  function _debounce(func, wait, immediate) {
    var _timeout;

    return function () {

      var context = this,
        args = arguments,
        later = function () {
          _timeout = null;
          if (!immediate) func.apply(context, args);
        };

      var callNow = immediate && !_timeout;

      clearTimeout(_timeout);

      _timeout = setTimeout(later, wait);

      if (callNow) {
        func.apply(context, args);
      }
    };
  }

  function _resetSearch() {
    $inputBox.val('');
    $resultsContainer.hide();
  }

  $.get(base_url + '/mkdocs/search_index.json', function (data) {
    var index = lunr(function () {
      this.field('title', {
        boost: 10
      });
      this.field('text');
      this.ref('location');
    });

    var documents = {};
    var doc;

    for (var i = 0; i < data.docs.length; i++) {
      doc = data.docs[i];
      doc.location = base_url + doc.location;
      index.add(doc);
      documents[doc.location] = doc;
    }

    function _search() {
      var query = $.trim($inputBox.val());
      $results.empty();

      console.log(query);

      if (query.length < _VALID_QUERY_LENGTH || query === '') {
        $resultsContainer.hide();
        return;
      }

      var results = index.search(query);
      var resultsHTML = '';
      var resultLength = results.length || null;

      if (query.length >= _VALID_QUERY_LENGTH && results.length === 0) {
        resultLength = 0;
      }

      if (resultLength !== null) {
        $resultsContainer.show();
        $searchContentBody.html('<strong><i class="fa fa-file-o"></i> ' + resultLength + '</strong> results found for <strong>' + query + '</strong>');
      }

      console.log("end of search function")
    }

    searchInput.addEventListener('keyup', _debounce(_search, 300));
  });

});