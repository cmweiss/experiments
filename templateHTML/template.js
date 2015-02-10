function loadTemplate() {
	var body = document.getElementsByTagName("body")[0];
	var midsec = document.createElement("div");
	var content = document.createElement("div");
	var header = document.createElement("div");
	var footer = document.createElement("div");
	var leftSide = document.createElement("div");
	var rightSide = document.createElement("div")
		bodyContents = document.createRange();
	
	header.setAttribute("id","header");
	midsec.setAttribute("id","midsection");
	content.setAttribute("id","content");
	footer.setAttribute("id","footer");
	leftSide.setAttribute("id","leftside");
	rightSide.setAttribute("id","rightside");
	
	bodyContents.setStart(body.firstChild, 0);
	bodyContents.setEndBefore(body.getElementsByTagName('script')[0]);
	
	content.appendChild(bodyContents.extractContents());
	body.appendChild(content);
	body.insertBefore(midsec,content);
	body.insertBefore(header,midsec);
	midsec.appendChild(leftSide);
	midsec.appendChild(content);
	midsec.appendChild(rightSide);
	body.appendChild(footer);
	ajaxFunction("header.htm","header");
	ajaxFunction("leftnav.htm","leftside");
	ajaxFunction("rightnav.htm","rightside");
	ajaxFunction("footer.htm","footer");
}
