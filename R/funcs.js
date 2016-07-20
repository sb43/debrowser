shinyjs.getNames = function(){
    var count = document.getElementsByClassName('tick').length; 
    var start = 0;  
    while(document.getElementsByClassName('tick')[start].getElementsByTagName('line')[0].getAttribute('x2') == 0){
    start += 1;
    }  
    var out = ''; 
    for (i = start; i < count; i++) 
    { 
    if('opacity: 1;' == document.getElementsByClassName('tick')[i].getAttribute('style')){ 
    out += document.getElementsByClassName('tick')[i].getElementsByTagName('text')[0].innerHTML + ',';
    }
    } 
    //document.getElementById('genenames').innerHTML = out;
    Shiny.onInputChange('genenames', out);
}