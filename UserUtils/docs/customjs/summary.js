var summary_text = $('.summary').html() ;
summary_text = summary_text.replace(/\|/g,'');
summary_text = "<span>Contents</span>" + summary_text;
$('.summary').html(summary_text);