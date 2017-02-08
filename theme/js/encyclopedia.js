$(document).ready(()=>{
    $('.header-text-link').attr('href','/Encyclopedia/');
})
$(document).on('click', '.markdown-link', e => {
    var el = $('.encyclopedia-content');
    var link = e.target.id;
    el.empty();
    fetch(link)
        .then(response => response.text())
        .then(text => {
            var converter = new showdown.Converter();
            var html = converter.makeHtml(text);
            el.append(html);
        })
        return false;
})
$(document).on('click', '.term-link, #view-all-link', e => {
    console.log($(e.target).attr("class"));
    var el = $('.encyclopedia-content');
    var targetId = $(e.target).attr('id');
    var hasTargetClass = $(e.target).hasClass('term-link');
    el.empty();

    fetch('/Encyclopedia/pages.json')
        .then(response => response.json())
        .then(json => {
            _(json)
                // .pickBy((value, key) => ((key[0]===targetId[0]) || (targetId.slice(0,3) === 'NUM' && !NaN(key[0]))) && targetClass === 'link-button' || targetId === 'view-all-link')
                .pickBy((value, key) => {
                    const isPickAChar = key[0] === targetId[0];
                    const wasPickANum = targetId.slice(0,3) === 'NUM';
                    const termStartsWithNumber = !isNaN(key[0]);
                    const wasTermLinkClicked = hasTargetClass;
                    const wasViewAllLinkClicked = targetId === 'view-all-link';
                    return (((isPickAChar || (wasPickANum && termStartsWithNumber)) && wasTermLinkClicked) || wasViewAllLinkClicked);
                })
                .forEach((value, key) => {
                    var keyPrint = key.replaceAll('_',' ');
                    el.append(`<div>
                        <a class='markdown-link' href='#' id='/Encyclopedia/Pages/Text/${value}'>${keyPrint}</a>
                    </div>`);
                })

            
            if(hasTargetClass){
                var first =(targetId[0]).charCodeAt(0) - 1;//previous term
                var last = (targetId[0]).charCodeAt(0) + 1;//next term
                if($('.adj-term-link-container').length == 0){
                    $('.encyclopedia-container').append(`
                        <div class = "adj-term-link-container" style="display: flex; flex-direction: row""></div>
                        `);//add container to main one to store the pointers
                }
                el2 = $('.adj-term-link-container');
                el2.empty();
                if(first > 64){
                    el2.append(`<div style="margin-right: auto">
                        <a class="term-link" href="#" id="${String.fromCharCode(first)}-button">
                            <i class="fa fa-chevron-left"></i>
                            Previous: \"${String.fromCharCode(first)}\" Terms
                        </a>
                    </div>`)
                }
                if(first === 64){//makes previous point to number
                    el2.append(`<div style="margin-right: auto">
                    
                    <a class="term-link" href="#" id="NUM-button">
                        <i class="fa fa-chevron-left"></i> 
                        Previous: \"#\" Terms
                    </a>
                    </div>`)
                }
                if(last < 91){
                    el2.append(`<div style="margin-left: auto">
                        <a class="term-link" href="#"id="${String.fromCharCode(last)}-button">
                            Next: \"${String.fromCharCode(last)}\" Terms
                            <i class="fa fa-chevron-right"></i>
                        </a>
                    </div>`)
                }
                if(last === 91){//makes next point to number
                        el2.append(`<div style="margin-left: auto">
                        
                        <a class="term-link" href="#"id="NUM-button">
                            Next: \"#\" Terms
                            <i class="fa fa-chevron-right"></i>
                        </a>
                        
                    </div>`)
                }
            }
            if(el.is(':empty')){
                el.append('<div>No documents match your query</div>');
            }
        })
}); 

String.prototype.replaceAll = function(search, replacement) {
    var target = this;
    return target.replace(new RegExp(search, 'g'), replacement);
};

