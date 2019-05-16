$(function () {
    var $inputBox = $('#hp-search__input');
    var searchInput = document.getElementById('hp-search__input');
    var $submit = $('#hp-search__submit');
    var $resultsContainer = $('#hp-search__results-container');
    var $results = $('#hp-search__results');
    var _VALID_QUERY_LENGTH = 3;
    var _isSearchActive = false;

    searchInput.addEventListener('keyup', _debounce(_search, 300));

    function _search() {
        var query = $.trim($inputBox.val());
        $results.empty();

        console.log(query);

        if (query.length < _VALID_QUERY_LENGTH || query === '') {
            $resultsContainer.hide();
            return;
        }
    }

    function _resetSearch() {
        $inputBox.val('');
        $resultsContainer.hide();
    }

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
});