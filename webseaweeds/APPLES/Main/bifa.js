var document = null;
var root = null;
//var ns = 'http://www.w3.org/2000/';
//var xlinkns = 'http://www.w3.org/1999/xlink';
var toolTip = null;
var trueCoords = null;
var tipBox = null;
var tipText = null;
var tipTitle = null;
var tipDesc = null;
var tipFactors = null;
var tipIupac = null;
var tipSeq = null;

var lastElement = null;
var titleText = '';
var titleDesc = '';
var titleFactors = '';
var titleIupac = '';
var titleSeq = '';


function Init(evt)
{
	document = evt.target.ownerDocument;
	root = document.documentElement;
	trueCoords = root.createSVGPoint();
	
	toolTip = document.getElementById('ToolTip');
	tipBox = document.getElementById('tipbox');
	tipText = document.getElementById('tipText');
	tipTitle = document.getElementById('tipTitle');
	tipDesc = document.getElementById('tipDesc');
	tipFactors = document.getElementById('tipFactors');
	tipIupac = document.getElementById('tipIupac');
	tipSeq = document.getElementById('tipSeq');
};


function GetTrueCoords(evt)
{
	// find the current zoom level and pan setting, and adjust the reported
	//    mouse position accordingly
	var newScale = root.currentScale;
	var translation = root.currentTranslate;
	trueCoords.x = (evt.clientX - translation.x)/newScale;
	trueCoords.y = (evt.clientY - translation.y)/newScale;
};




function ToggleDisplay(el)
{
	var current = el.getAttributeNS(null, 'display');
	if ('none' == current)
	{
		el.setAttributeNS(null, 'display', 'inline');
		return true;
	}
	else
	{
		el.setAttributeNS(null, 'display', 'none');
		return false;
	}
};



function FactorMouseOver(factor_id, show)
{
	var index = null;
	var factor = document.getElementById(factor_id);
	if (null != factor)
	{
		var children = factor.getElementsByTagName('bifa:hit');
		for(var i=0; i < children.length; i++)
		{
			idx = children.item(i).getAttributeNS('', 'index');
			var bg_id = 'hit_' + idx + '_bg';
			var bg = document.getElementById(bg_id);
			if (null != bg)
			{
				if (show)
				{
					bg.setAttributeNS('', 'display', 'inline');
				}
				else
				{
					bg.setAttributeNS('', 'display', 'none');
				}
			}
		}
	}
}



function NotesMouseOver( show )
{
	var index = null;
	var panel = document.getElementById( 'notes_panel' );
	if ( null != panel )
	{
		if (show)
		{
			var back = document.getElementById( 'notes_background' );
			var outline = panel.getBBox();
			back.setAttributeNS(null, 'width', Number(outline.width) + 2);
			back.setAttributeNS(null, 'height', Number(outline.height) + 2);
			back.setAttributeNS(null, 'x', Number(outline.x) - 1);
			back.setAttributeNS(null, 'y', Number(outline.y) - 1);
			panel.setAttributeNS( '', 'display', 'inline' );
		}
		else
		{
			//alert( 'hello' );
			panel.setAttributeNS( '', 'display', 'none' );
		}
	}
}



function OnClick(evt)
{
	try
	{
		var id = null;
		var element = evt.target;
		do
		{
			id = element.getAttributeNS(null, 'id');
			if ('' != id)
			{
				break;
			}
			element = element.parentNode;
		}
		while (null != element)
			
			//alert(id);
		
		var info = document.getElementById(id + '_info');
		if (null != info)
		{
			var parent = element.parentNode;
			parent.removeChild(element);
			parent.appendChild(element);
			
			var back = document.getElementById(id + '_info_back');
			var underlines = document.getElementById(id + '_underlines');
			if (ToggleDisplay(info))
			{
				underlines.setAttributeNS(null, 'display', 'inline');
				var outline = info.getBBox();
				back.setAttributeNS(null, 'width', Number(outline.width) + 2);
				back.setAttributeNS(null, 'height', Number(outline.height) + 2);
				back.setAttributeNS(null, 'x', Number(outline.x) - 1);
				back.setAttributeNS(null, 'y', Number(outline.y) - 1);
			}
			else
			{
				underlines.setAttributeNS(null, 'display', 'none');
				back.setAttributeNS(null, 'width', 0);
				back.setAttributeNS(null, 'height', 0);
				back.setAttributeNS(null, 'x', 0);
				back.setAttributeNS(null, 'y', 0);
			}
		}
		
	}
	catch (er)
	{
	}
};



function ShowTooltip(evt, turnOn)
{
	try
	{
		if (!evt || !turnOn)
		{
			toolTip.setAttributeNS(null, 'display', 'none');
		}
		else
		{
			var tipScale = 1/root.currentScale;
			var textWidth = 0;
			var tspanWidth = 0;
			var boxHeight = 20;
			
			tipBox.setAttributeNS(null, 'transform', 'scale(' + tipScale + ',' + tipScale + ')' );
			tipText.setAttributeNS(null, 'transform', 'scale(' + tipScale + ',' + tipScale + ')' );
			
			var targetElement = evt.target;
			
			if ( lastElement != targetElement )
			{
				var targetTitle = targetElement.getElementsByTagName('title').item(0);
				titleText = null;
				if ( targetTitle )
				{
					titleText = targetTitle.firstChild.nodeValue;
				}
				tipTitle.firstChild.nodeValue = titleText;
				
				var targetDesc = targetElement.getElementsByTagName('desc').item(0);
				titleDesc = null;
				if ( targetDesc )
				{
					titleDesc = targetDesc.firstChild.nodeValue;
				}
				tipDesc.firstChild.nodeValue = titleDesc;
				
				var targetFactors = targetElement.getElementsByTagName('factors').item(0);
				titleFactors = null;
				if ( targetFactors )
				{
					titleFactors = targetFactors.firstChild.nodeValue;
				}
				tipFactors.firstChild.nodeValue = titleFactors;
				
				var targetIupac = targetElement.getElementsByTagName('iupac').item(0);
				titleIupac = null;
				if ( targetIupac )
				{
					titleIupac = targetIupac.firstChild.nodeValue;
				}
				tipIupac.firstChild.nodeValue = titleIupac;
				
				var targetSeq = targetElement.getElementsByTagName('seq').item(0);
				titleSeq = null;
				if ( targetSeq )
				{
					titleSeq = targetSeq.firstChild.nodeValue;
				}
				tipSeq.firstChild.nodeValue = titleSeq;
				
			}
			
			if ( titleText || titleDesc )
			{
				var xPos = trueCoords.x;
				var yPos = trueCoords.y;
				
				//return rectangle around object as Rect object
				var outline = tipText.getBBox();
				tipBox.setAttributeNS(null, 'width', Number(outline.width) + 10);
				tipBox.setAttributeNS(null, 'height', Number(outline.height) + 10);
				
				//toolTip.setAttributeNS(null, 'transform', 'translate(' + xPos + ',' + yPos + ')');
				toolTip.setAttributeNS(null, 'display', 'inline');
			}
		}
	}
	catch(er){}
};




