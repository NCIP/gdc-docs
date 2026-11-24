---
template: encyclopedia.html
---

<script src="https://code.jquery.com/jquery-3.7.1.min.js" integrity="sha256-/JqT3SQfawRcv/BIHPThkBvs0OEvtFFmqPF/lYI/Cxo=" crossorigin="anonymous"></script>
<script type="text/javascript" src="js/encyclopedia.js"></script>
<link rel="stylesheet" href="css/encyclopedia.css">
<script>
window.navigation.addEventListener("navigate", (event) => {
   updateDictView(event.destination.url); 
});
</script>
<script>
window.addEventListener("load", () => {
    updateDictView(window.location.href);
});
</script>

# GDC Encyclopedia
