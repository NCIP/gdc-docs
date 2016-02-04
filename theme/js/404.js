$(function() {
  'use strict';

  function Redirector(redirectJSONURL, timeOutMS) {

    var _self = this,
        _currentPath = (location.pathname + location.hash).toLowerCase(),
        _redirectMap = null,
        _timeOutMS = timeOutMS || 5000,
        _redirectJSONConfigURL = redirectJSONURL || '/config/redirects.json';

    function _performRedirect(defer) {
      var redirects = _redirectMap.redirects;

      for (var i = 0; i < redirects.length; i++) {
        var redirectRule = redirects[i],
            fromURLRedirect = redirectRule.from.toLowerCase();

        // Note this is a case insensitive (exact) string match trailing slashes may break the comparison
        if (_currentPath.indexOf(fromURLRedirect) >= 0 && typeof redirectRule.to === 'string') {
          var targetURL = _redirectMap.defaultBaseToURL + redirectRule.to;

          defer.resolve(targetURL);

          setTimeout(function(){
            window.location.href = _redirectMap.defaultBaseToURL + redirectRule.to;
          }, _timeOutMS);

          return;
        }
      }

      defer.reject(null);
    }

    function _loadRedirectConfig(defer) {

      var loadConfigPromise = $.Deferred();;

      $.get(_redirectJSONConfigURL)
        .done(function(configRedirectData) {
          _redirectMap = configRedirectData;
        })
        .fail(function(error) {
          defer.reject(error);
        })
        .always(function() {
          loadConfigPromise.resolve();
        });

      return loadConfigPromise;
    }

    function _redirect() {

      var redirectPromise = $.Deferred();

      if (_redirectMap) {
        _performRedirect(redirectPromise);
      }
      else {
        _loadRedirectConfig()
          .then(function() {
            _performRedirect(redirectPromise);
          });
      }

      return redirectPromise;
    }


    _self.redirect = _redirect;
  }

  var TIMEOUT_IN_SECONDS = 5,
      redirector = new Redirector(null, (TIMEOUT_IN_SECONDS * 1000)),
      pageNotFoundContainer = $('#page-not-found-container'),
      pageFoundContainer = $('#redirect-page-found-container'),
      redirectLinkEl = $('#redirect-link'),
      redirectCountdownEl = $('#redirect-timer'),
      loadingEl = $('#redirect-lookup-container');

  function renderTimer(timeoutInSeconds) {

    var text = timeoutInSeconds;

    if (timeoutInSeconds > 0) {
      redirectCountdownEl.stop().animate({opacity: 1}, 100);
    }
    else {
      redirectCountdownEl.stop().css({opacity: 1}, 100);
    }


    if (timeoutInSeconds === 0) {
      text = '0 (redirecting)';
    }

    redirectCountdownEl.text(text);

    if (timeoutInSeconds > 0) {
      redirectCountdownEl.animate({opacity: 0}, 900);
    }

    timeoutInSeconds--;

    if (timeoutInSeconds >= 0) {
      setTimeout(function() { renderTimer(timeoutInSeconds); }, 1000);
    }
  }



  redirector
    .redirect()
    .then(function(destinationURL) {
      loadingEl.hide();
      pageNotFoundContainer.hide();
      redirectLinkEl.attr('href', destinationURL).text(destinationURL);
      pageFoundContainer.show();
      renderTimer(TIMEOUT_IN_SECONDS);
    })
    .fail(function() {
      loadingEl.hide();
      pageNotFoundContainer.show();
    });

});