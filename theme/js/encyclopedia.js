$(function () {
    url = window.location.href;
    if(url.endsWith('/')){
        $('.doc-list').addClass('hidden');
        $('.encyclopedia.content').removeClass('hidden');
    }
    if(url.endsWith('/#')){
        $('.encyclopedia-content').addClass('hidden');
        $('.doc-list').removeClass('hidden');
    }
    if((url.slice(-1).charCodeAt(0) >= 65 && url.slice(-1).charCodeAt(0) <= 90) || url.slice(-5) === '/#num'){
        $('.encyclopedia-content').addClass('hidden');
        if(url.slice(-5) == '/#num'){
            $('.doc-list').not('#num-list').addClass('hidden');
            $('#num-list').removeClass('hidden');
        }
        else{
            $('#' + url.slice(-1) + '-list').removeClass('hidden');
            $('.doc-list').not('#' + url.slice(-1) + '-list').addClass('hidden');
            $('#num-list').addClass('hidden');
        }

    }
});

