# app.py

import streamlit as st
from data import research_summaries
from collections import defaultdict

# -------------------------
# Set Page Configuration
# -------------------------
st.set_page_config(
    page_title="📚 Bioinformatics Learning App",
    layout="wide",
    initial_sidebar_state="expanded",
)

# -------------------------
# Custom CSS for Styling (Optional)
# -------------------------
def local_css(file_name):
    try:
        with open(file_name) as f:
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
    except FileNotFoundError:
        st.warning("CSS file not found. Skipping custom styling.")

# Uncomment the following line if you have a CSS file
# local_css("assets/styles.css")

# -------------------------
# Initialize Session State
# -------------------------
if 'current_index' not in st.session_state:
    st.session_state.current_index = 0

if 'likes' not in st.session_state:
    st.session_state.likes = defaultdict(bool)

if 'comments' not in st.session_state:
    st.session_state.comments = defaultdict(list)

if 'selected_category' not in st.session_state:
    st.session_state.selected_category = "All"

# -------------------------
# Sidebar: Category Selection
# -------------------------
st.sidebar.title("📖 Categories")

# Extract unique categories
categories = sorted(list(set([item['category'] for item in research_summaries])))

# Add 'All' option
categories.insert(0, "All")

# Single select for category
selected_category = st.sidebar.selectbox(
    "Select a Category",
    options=categories,
    index=0,
    key='category_select'
)

st.session_state.selected_category = selected_category

# -------------------------
# Main Title
# -------------------------
st.title("📚 Bioinformatics Learning App")

# -------------------------
# Search Functionality
# -------------------------
search_query = st.text_input(
    "🔍 Search Topics",
    value="",
    placeholder="Enter keywords to search...",
    help="Search within titles, summaries, and tags."
)

# -------------------------
# Filter Summaries Based on Category and Search
# -------------------------
def filter_summaries(summaries, category, query):
    filtered = summaries.copy()
    if category != "All":
        filtered = [s for s in filtered if s['category'] == category]
    if query:
        query = query.lower()
        filtered = [
            s for s in filtered
            if query in s['title'].lower()
            or query in s['content'].lower()
            or any(query in tag.lower() for tag in s['tags'])
        ]
    return filtered

filtered_summaries = filter_summaries(research_summaries, st.session_state.selected_category, search_query)

# -------------------------
# Handle No Results
# -------------------------
if not filtered_summaries:
    st.warning(f"No results found for '{search_query}' in category '{selected_category}'.")
    st.stop()

# -------------------------
# Pagination: Current Summary
# -------------------------
current_summary = filtered_summaries[st.session_state.current_index]

# -------------------------
# Display Summary and Visualizations
# -------------------------
with st.container():
    st.subheader(current_summary['title'])
    st.markdown(current_summary['content'])
    st.markdown(f"**Category:** {current_summary['category']}")
    st.markdown(f"**Tags:** {', '.join(current_summary['tags'])}")

    # -------------------------
    # Action Buttons
    # -------------------------
    col1, col2, col3 = st.columns(3)

    with col1:
        like_button = st.button(
            "❤️ Like",
            key=f"like_{current_summary['id']}",
            on_click=lambda: toggle_like(current_summary['id']),
        )
        if st.session_state.likes[current_summary['id']]:
            st.markdown("**You liked this!**")

    with col2:
        comment_input = st.text_input(
            "💬 Add a comment...",
            key=f"comment_input_{current_summary['id']}",
            placeholder="Type your comment here..."
        )
        if st.button("Submit Comment", key=f"submit_comment_{current_summary['id']}"):
            if comment_input:
                st.session_state.comments[current_summary['id']].append(comment_input)
                st.success("Comment submitted!")
            else:
                st.warning("Please enter a comment before submitting.")

    with col3:
        share_button = st.button("🔗 Share", key=f"share_{current_summary['id']}")
        if share_button:
            st.info("Share functionality is not implemented yet.")

    # -------------------------
    # Display Comments
    # -------------------------
    if st.session_state.comments[current_summary['id']]:
        st.markdown("### 💬 Comments")
        for idx, comment in enumerate(st.session_state.comments[current_summary['id']], 1):
            st.write(f"{idx}. {comment}")

# -------------------------
# Navigation Buttons
# -------------------------
st.markdown("---")
nav_col1, nav_col2, nav_col3 = st.columns([1, 2, 1])

with nav_col1:
    if st.button("⬅️ Previous", key="prev_button"):
        if st.session_state.current_index > 0:
            st.session_state.current_index -= 1
        else:
            st.warning("You're at the first topic.")

with nav_col3:
    if st.button("Next ➡️", key="next_button"):
        if st.session_state.current_index < len(filtered_summaries) - 1:
            st.session_state.current_index += 1
        else:
            st.warning("You've reached the last topic.")

# -------------------------
# Pagination Indicator
# -------------------------
st.markdown(
    f"**Topic {st.session_state.current_index + 1} of {len(filtered_summaries)}**"
)

# -------------------------
# Toggle Like Function
# -------------------------
def toggle_like(summary_id):
    st.session_state.likes[summary_id] = not st.session_state.likes[summary_id]
